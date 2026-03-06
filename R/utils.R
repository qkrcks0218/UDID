## ============================================================
## SHARED UTILITIES (internal)
## ============================================================

#' @keywords internal
safe_ratio <- function(num, den, clip = 1e-6) {
  num / pmax(den, clip)
}

#' @keywords internal
Truncate.Function <- function(X.Input, LB, UB) {
  if (is.null(dim(X.Input))) {
    apply(cbind(X.Input, LB, UB), 1, median)
  } else {
    LL  <- list(as.matrix(X.Input),
                matrix(LB, nrow(X.Input), ncol(X.Input)),
                matrix(UB, nrow(X.Input), ncol(X.Input)))
    ULL <- array(unlist(LL), c(nrow(X.Input), ncol(X.Input), 3))
    out <- apply(ULL, c(1, 2), median)
    colnames(out) <- colnames(X.Input)
    out
  }
}

#' @keywords internal
Mboot <- function(V, N, NumBoot = N) {
  BootMat1 <- matrix(rnorm(N * NumBoot) + 1, N, NumBoot)
  apply(BootMat1 * matrix(V, N, NumBoot), 2, mean)
}

#' @keywords internal
EIF_Continuous <- function(Y1, Y0, A, OR, bA1, mu1, DR, efficient = TRUE) {
  Eff.Ind <- ifelse(efficient, 1, 0)
  (1 - A) * bA1 * OR * (Y1 - mu1) +
    A * mu1 +
    Eff.Ind * (2 * A - 1) * DR * (Y0 - mu1)
}

## ---- Sensitivity scaling: alpha_0 -> alpha_1 ----------------
##
##  New approach (Park & Tchetgen Tchetgen, revised):
##  1. Compute mu(x) = E[Y1 | A=0, X=x] under Gamma=1 (baseline).
##  2. For Gamma > 1 define directional scales:
##
##     scale_UB(y, x) = Gamma^{1/2}   if y >  mu(x)
##                       Gamma^{-1/2}  if y <= mu(x)
##
##     scale_LB(y, x) = Gamma^{-1/2}  if y >  mu(x)
##                       Gamma^{1/2}   if y <= mu(x)
##
##  3. alpha_{1,UB} = alpha_0 * scale_UB
##     alpha_{1,LB} = alpha_0 * scale_LB
##
##  Works element-wise on vectors/matrices; y and mu must be
##  conformable (same length, or y is a matrix and mu a vector
##  recycled by row).

#' @keywords internal
sens_scale_UB <- function(y, mu, Gamma) {
  ## y and mu must be conformable; returns same shape as y
  g_up   <- Gamma^(0.5)
  g_down <- Gamma^(-0.5)
  ifelse(y > mu, g_up, g_down)
}

#' @keywords internal
sens_scale_LB <- function(y, mu, Gamma) {
  g_up   <- Gamma^(0.5)
  g_down <- Gamma^(-0.5)
  ifelse(y > mu, g_down, g_up)
}

## ---- Box-Cox helpers (parametric only) ----------------------

#' @keywords internal
bc_find_shift <- function(y) {
  if (any(y <= 0)) -min(y) + 1e-3 else 0
}

#' @keywords internal
bc_transform <- function(y, lambda) {
  if (abs(lambda) < 1e-10) log(y) else (y^lambda - 1) / lambda
}

#' @keywords internal
bc_inverse <- function(t, lambda) {
  if (abs(lambda) < 1e-10) exp(t) else (lambda * t + 1)^(1 / lambda)
}

## ---- Logistic density ratio classifier (parametric only) ----
##
##  Estimates  r(y,x) = f_num(y,x) / f_den(y,x)  via:
##    r(y,x) = [p1/(1-p1)] * [N_den / N_num]

#' @keywords internal
fit_density_ratio <- function(y_num, X_num, y_den, X_den) {
  N_num       <- length(y_num)
  N_den       <- length(y_den)
  ratio_prior <- N_den / max(N_num, 1)
  xnms        <- paste0("X", seq_len(ncol(X_num)))
  
  df  <- data.frame(label = c(rep(1, N_num), rep(0, N_den)),
                    y     = c(y_num, y_den),
                    rbind(X_num, X_den))
  colnames(df)[-(1:2)] <- xnms
  
  list(fit         = stats::glm(label ~ ., data = df, family = stats::binomial()),
       ratio_prior = ratio_prior,
       xnames      = xnms)
}

#' @keywords internal
predict_density_ratio <- function(dr_obj, y_eval, X_eval) {
  df_eval <- data.frame(y = y_eval, X_eval)
  colnames(df_eval)[-1] <- dr_obj$xnames
  p1  <- pmin(pmax(stats::predict(dr_obj$fit, newdata = df_eval, type = "response"), 1e-6), 1 - 1e-6)
  pmin(pmax((p1 / (1 - p1)) * dr_obj$ratio_prior, 1e-4), 1e4)
}

#' @keywords internal
fit_density_ratio_NoX <- function(y_num, y_den) {
  N_num       <- length(y_num)
  N_den       <- length(y_den)
  ratio_prior <- N_den / max(N_num, 1)
  
  df  <- data.frame(label = c(rep(1, N_num), rep(0, N_den)),
                    y     = c(y_num, y_den))
  
  list(fit         = stats::glm(label ~ y, data = df, family = stats::binomial()),
       ratio_prior = ratio_prior)
}

#' @keywords internal
predict_density_ratio_NoX <- function(dr_obj, y_eval) {
  df_eval <- data.frame(y = y_eval)
  p1  <- pmin(pmax(stats::predict(dr_obj$fit, newdata = df_eval, type = "response"), 1e-6), 1 - 1e-6)
  pmin(pmax((p1 / (1 - p1)) * dr_obj$ratio_prior, 1e-4), 1e4)
}

## ---- Compile Rcpp kernels -----------------------------------
#sourceCpp("udid_kliep.cpp")   # compiles gaussian_kernel_ratio()

## ---- KLIEP kernel density ratio (nonparametric only) --------
## Replaces the nested-apply R loop in densratio:::compute_kernel_Gaussian
## with a parallelised C++ implementation. Drop-in compatible.

#' @keywords internal
compute_density_ratio_Kernel <- function(x, FIT.LIST) {
  if (is.vector(x)) x <- matrix(x)
  sigma <- FIT.LIST$kernel_info$sigma
  if (length(sigma) > 1) {
    ## Per-dimension bandwidth
    gaussian_kernel_ratio_bw(
      x              = as.matrix(x),
      centers        = as.matrix(FIT.LIST$kernel_info$centers),
      sigma_vec      = as.numeric(sigma),
      kernel_weights = as.numeric(FIT.LIST$kernel_weights)
    )
  } else {
    ## Scalar bandwidth (original path)
    gaussian_kernel_ratio(
      x              = as.matrix(x),
      centers        = as.matrix(FIT.LIST$kernel_info$centers),
      sigma          = sigma,
      kernel_weights = as.numeric(FIT.LIST$kernel_weights)
    )
  }
}

## ---- KLIEP with per-dimension bandwidths (sigma_Y, sigma_X) ----
##
## Replaces densratio::densratio(..., method = "KLIEP") when different
## bandwidths are needed for the Y columns and the X columns.
##
## Usage:
##   fit <- KLIEP_bw(x1, x2, n_y_cols = 1)
##   r   <- compute_density_ratio_Kernel(eval_pts, fit)
##
## The sigma search is a coordinate-descent grid over (sigma_Y, sigma_X),
## using cross-validated KLIEP log-likelihood, identical to the original
## densratio:::KLIEP_search_sigma but in two dimensions.

#' @keywords internal
KLIEP_bw <- function(x1, x2, n_y_cols = 1, sigma = "median",
                     kernel_num = 100, fold = 5, verbose = FALSE) {
  if (is.vector(x1)) x1 <- matrix(x1)
  if (is.vector(x2)) x2 <- matrix(x2)
  x1 <- as.matrix(x1); x2 <- as.matrix(x2)
  
  nx1 <- nrow(x1)
  p   <- ncol(x1)
  kernel_num <- min(kernel_num, nx1)
  centers <- x1[sample(nx1, size = kernel_num), , drop = FALSE]
  
  ## --- Determine sigma_vec and fit weights ---
  sigma_init <- .median_heuristic_bw(x1, x2, n_y_cols)
  
  if (identical(sigma, "median")) {
    ## Use median heuristic directly, just optimise weights
    sigma_vec <- sigma_init
    phi_x1 <- gaussian_kernel_matrix_bw(x1, centers, sigma_vec)
    phi_x2 <- gaussian_kernel_matrix_bw(x2, centers, sigma_vec)
    kernel_weights <- as.vector(densratio:::KLIEP_optimize_alpha(phi_x1, phi_x2))
  } else {
    ## "auto": full CV search + weight optimisation, all in C++
    res <- kliep_bw_cpp(x1, x2, centers, sigma_init*10, n_y_cols, fold)
    sigma_vec      <- as.numeric(res$sigma)
    kernel_weights <- as.numeric(res$kernel_weights)
    centers        <- as.matrix(res$centers)
  }
  
  list(
    kernel_weights = kernel_weights,
    kernel_info    = list(kernel     = "Gaussian",
                          kernel_num = kernel_num,
                          sigma      = sigma_vec,
                          centers    = centers),
    fold           = fold
  )
}

## Median heuristic for per-dimension bandwidths:
##   sigma_Y = median of pairwise distances among Y columns
##   sigma_X = median of pairwise distances among X columns
#' @keywords internal
.median_heuristic_bw <- function(x1, x2, n_y_cols) {
  p     <- ncol(x1)
  idx_Y <- seq_len(n_y_cols)
  idx_X <- if (n_y_cols < p) (n_y_cols + 1):p else integer(0)
  
  idx_sub <- 1:nrow(x1)
  
  sigma_Y <- stats::median(stats::dist(x1[idx_sub, idx_Y, drop = FALSE]))
  if (sigma_Y < 1e-6) sigma_Y <- 1.0
  
  if (length(idx_X) > 0) {
    sigma_X <- stats::median(stats::dist(x1[idx_sub, idx_X, drop = FALSE]))
    if (sigma_X < 1e-6) sigma_X <- 1.0
  } else {
    sigma_X <- 1.0
  }
  
  sv <- numeric(p)
  sv[idx_Y] <- sigma_Y
  if (length(idx_X) > 0) sv[idx_X] <- sigma_X
  sv
}
 