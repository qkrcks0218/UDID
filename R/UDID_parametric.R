#' Parametric UDID Estimator (Experimental/Beta Version)
#'
#' Computes the average treatment effect on the treated (ATT)
#' under the universal difference-in-differences (UDID) framework
#' (Park and Tchetgen Tchetgen, 2026+). Nuisance functions are
#' estimated parametrically based on parametric models. Sensitivity
#' analysis is supported for user-specified sensitivity parameter values.
#'
#' @param Y0 Numeric vector of pre-treatment outcomes.
#' @param Y1 Numeric vector of post-treatment outcomes.
#' @param A Binary treatment indicator (0/1).
#' @param X Optional covariate matrix. If \code{NULL}, a no-covariate version is used.
#' @param type Outcome type: \code{"continuous"}, \code{"binary"}, or \code{"poisson"}.
#' @param log_Gamma_seq Numeric scalar or vector of log(Gamma) sensitivity values (default 0).
#'
#' @details
#' \code{UDID_Parametric} implements the parametric universal
#' difference-in-differences method. 
#' Nuisance functions are fitted once based on parametric models. 
#' The ATT estimate is computed for every value of the sensitivity parameter \eqn{\Gamma} 
#' at negligible additional cost.
#'
#' The key assumption is odds ratio equi-confounding (OREC), which states that
#' \eqn{\alpha_1(y,x) = \alpha_0(y,x)}, where
#' \deqn{
#'   \alpha_t(y,x) =
#'   \frac{f(Y_t^{(0)}=y \mid A=1, X=x)}{f(Y_t^{(0)}=y_R \mid A=1, X=x)}
#'   \frac{f(Y_t^{(0)}=y_R \mid A=0, X=x)}{f(Y_t^{(0)}=y \mid A=0, X=x)}
#' }
#' and \eqn{y_R} is a reference value. Note that \eqn{\alpha_t(y_R,X)=1}.
#'
#' When the outcome is \strong{continuous}, the nuisance functions are estimated
#' as follows:
#' \itemize{
#'   \item \eqn{f(Y_1 \mid A=0, X)} is modeled via a Box-Cox transformation
#'   (Box and Cox, 1964) followed by a Gaussian linear model.
#'   \item \eqn{f(Y_0 \mid A=1, X)/f(Y_0 \mid A=0, X)} is estimated via a
#'   logistic regression density ratio classifier, which classifies numerator
#'   versus denominator samples and converts predicted probabilities to density
#'   ratios.
#'   \item \eqn{f(A \mid X)} is estimated via logistic regression.
#' }
#'
#' When the outcome is \strong{binary}, all conditional distributions
#' (\eqn{f(A \mid X)}, \eqn{f(Y_0 \mid A, X)}, and \eqn{f(Y_1 \mid A=0, X)})
#' are estimated via logistic regression (GLM with binomial family).
#'
#' When the outcome is \strong{Poisson} (\code{type = "poisson"}), the outcome
#' regressions \eqn{f(Y_0 \mid A=0, X)} and \eqn{f(Y_1 \mid A=0, X)} are
#' estimated via Poisson GLM, \eqn{f(A \mid X)} via logistic regression, and
#' the odds ratio \eqn{\alpha_0} via a logistic density ratio classifier.
#'
#' Given these nuisance estimates, the ATT estimate is obtained via the efficient influence function.
#'
#' For sensitivity analysis, given the sensitivity parameter \eqn{\Gamma \geq 1},
#' we allow
#' \deqn{
#'   \alpha_1(y,x) \in
#'   \left[\Gamma^{-1} \cdot \alpha_0(y,x),\; \Gamma \cdot \alpha_0(y,x)\right].
#' }
#' The sensitivity bounds are computed using the same two-step approach
#' described in \code{\link{UDID_Nonparametric}}.
#'
#' @return A named list with the following components:
#'   \describe{
#'     \item{\code{Effect}}{A single-row data frame reporting the estimated
#'       average treatment effect on the treated (ATT) under the UDID method
#'       at \code{log_Gamma = 0} and its asymptotic standard error (\code{SE}).}
#'     \item{\code{Sensitivity}}{A data frame with one row per entry of
#'       \code{log_Gamma_seq}, reporting the lower and upper sensitivity bounds
#'       on the ATT (\code{ATT_LB}, \code{ATT_UB}) along with their asymptotic
#'       standard errors (\code{SE_LB}, \code{SE_UB}).}
#'   }
#'
#' @references
#' \itemize{
#'   \item Park, C., & Tchetgen Tchetgen, E. (2026+).
#'     A Universal Nonparametric Framework for Difference-in-Differences Analyses.
#'     \url{https://arxiv.org/abs/2212.13641}.
#'   \item Box, G. E. P., & Cox, D. R. (1964).
#'     An analysis of transformations.
#'     \emph{Journal of the Royal Statistical Society, Series B}, 26(2), 211--243.
#' }
#'
#' @seealso \code{\link{UDID_Nonparametric}} for the nonparametric version,
#'   \code{\link{UDID_Sensitivity_Bounds}} and \code{\link{UDID_Sensitivity_Plot}}
#'   for post-estimation sensitivity analysis tools.
#'
#' @export
UDID_Parametric <- function(Y0,
                            Y1,
                            A,
                            X             = NULL,
                            type          = "continuous",
                            log_Gamma_seq = 0) {
  
  ## Dispatch to no-covariate version if X is NULL
  if (is.null(X)) {
    return(UDID_Parametric_NoX(Y0 = Y0, Y1 = Y1, A = A,
                               type = type, log_Gamma_seq = log_Gamma_seq))
  }
  
  Gamma_seq <- exp(log_Gamma_seq)
  N      <- length(Y0)
  d      <- ncol(X)
  xnames <- paste0("X", sprintf("%05d", seq_len(d)))
  df_X   <- stats::setNames(as.data.frame(X), xnames)
  idx0   <- which(A == 0)
  idx1   <- which(A == 1)
  pr_A   <- mean(A)
  
  ## ----------------------------------------------------------
  ## Fit all nuisance models ONCE, independent of Gamma
  ## ----------------------------------------------------------
  
  if (type == "continuous") {
    
    ## f1(Y1 | A=0, X) via Box-Cox + Gaussian GLM
    Y1_A0     <- Y1[idx0]
    shift_Y1  <- bc_find_shift(Y1_A0)
    Y1_A0_pos <- Y1_A0 + shift_Y1
    
    bc_Y1     <- MASS::boxcox(Y1_A0_pos ~ as.matrix(df_X[idx0, ]),
                              lambda = seq(-2, 2, by = 0.05), plotit = FALSE)
    lambda_Y1 <- bc_Y1$x[which.max(bc_Y1$y)]
    
    T_Y1_A0   <- bc_transform(Y1_A0_pos, lambda_Y1)
    fit_f1    <- stats::lm(TY1 ~ ., data = cbind(TY1 = T_Y1_A0, df_X[idx0, , drop = FALSE]))
    sigma2_f1 <- stats::var(stats::residuals(fit_f1))
    
    ## pr(A=1 | X)
    fit_ps  <- stats::glm(A ~ ., data = cbind(A = A, df_X), family = stats::binomial())
    ps_eval <- pmin(pmax(stats::predict(fit_ps, newdata = df_X, type = "response"), 1e-6), 1 - 1e-6)
    odds_X  <- safe_ratio(ps_eval, 1 - ps_eval)
    
    ## Density ratio classifiers (Gamma-free)
    dr0 <- fit_density_ratio(Y0[idx1], X[idx1, , drop = FALSE],
                             Y0[idx0], X[idx0, , drop = FALSE])
    dr1 <- fit_density_ratio(Y1[idx0], X[idx0, , drop = FALSE],
                             Y0[idx0], X[idx0, , drop = FALSE])
    
    y_ref     <- stats::median(Y0[idx0])
    r0_Y0     <- predict_density_ratio(dr0, Y0,           X)
    r0_Y1     <- predict_density_ratio(dr0, Y1,           X)
    r0_ref    <- predict_density_ratio(dr0, rep(y_ref, N), X)
    
    ## alpha_0 (Gamma-free)
    alpha0_Y0 <- safe_ratio(r0_Y0, r0_ref)
    alpha0_Y1 <- safe_ratio(r0_Y1, r0_ref)
    beta0     <- r0_ref * odds_X
    r1_Y0     <- predict_density_ratio(dr1, Y0, X)
    
    ## Grid integration — alpha_0 on equally-spaced grid (Gamma-free)
    n_grid    <- 501
    mu1_bc    <- stats::predict(fit_f1, newdata = df_X)
    sig1_bc   <- sqrt(sigma2_f1)
    
    ## Equally-spaced grid in Box-Cox space, covering ±4 sd
    t_lo      <- min(mu1_bc) - 4 * sig1_bc
    t_hi      <- max(mu1_bc) + 4 * sig1_bc
    t_grid    <- seq(t_lo, t_hi, length.out = n_grid)
    dT        <- diff(t_grid)[1]
    
    ## Back-transform grid to original scale
    Y1_grid_raw   <- bc_inverse(t_grid, lambda_Y1) - shift_Y1
    valid_grid    <- is.finite(Y1_grid_raw)
    Y1_grid_raw[!valid_grid] <- y_ref
    
    ## Conditional density weights: f(t | X_i) * dT  for each (i, grid point)
    ## phi((t_j - mu_i) / sig) / sig * dT
    wt_mat <- outer(mu1_bc, t_grid,
                    FUN = function(mu, t) stats::dnorm(t, mean = mu, sd = sig1_bc)) * dT   # N x n_grid
    
    ## Evaluate alpha_0 on grid
    X_rep   <- X[rep(seq_len(N), times = n_grid), , drop = FALSE]
    y1_flat <- rep(Y1_grid_raw, each = N)
    
    alpha0_grid <- matrix(
      safe_ratio(predict_density_ratio(dr0, y1_flat,                X_rep),
                 predict_density_ratio(dr0, rep(y_ref, N * n_grid), X_rep)),
      N, n_grid)
    alpha0_grid[, !valid_grid] <- 0
    Y1_grid_mat <- matrix(Y1_grid_raw, N, n_grid, byrow = TRUE)
    Y1_grid_mat[, !valid_grid] <- 0
    
  } else if (type == "binary") {
    
    Y0r <- round(Y0); Y1r <- round(Y1); Ar <- round(A)
    df_all <- cbind(Y0 = Y0r, Y1 = Y1r, A = Ar, df_X)
    df_A0b <- df_all[Ar == 0, ]; df_A1b <- df_all[Ar == 1, ]
    
    fit_ps    <- stats::glm(A ~ ., data = df_all[, c("A", xnames)], family = stats::binomial())
    ps_eval   <- pmin(pmax(stats::predict(fit_ps, newdata = df_X, type = "response"), 1e-6), 1 - 1e-6)
    odds_X    <- safe_ratio(ps_eval, 1 - ps_eval)
    
    fit_Y0gA0 <- stats::glm(Y0 ~ ., data = df_A0b[, c("Y0", xnames)], family = stats::binomial())
    fit_Y0gA1 <- stats::glm(Y0 ~ ., data = df_A1b[, c("Y0", xnames)], family = stats::binomial())
    fit_Y1gA0 <- stats::glm(Y1 ~ ., data = df_A0b[, c("Y1", xnames)], family = stats::binomial())
    
    p_Y0_1_A0 <- stats::predict(fit_Y0gA0, newdata = df_X, type = "response")
    p_Y0_1_A1 <- stats::predict(fit_Y0gA1, newdata = df_X, type = "response")
    p_Y1_1_A0 <- stats::predict(fit_Y1gA0, newdata = df_X, type = "response")
    p_Y0_0_A0 <- 1 - p_Y0_1_A0;  p_Y0_0_A1 <- 1 - p_Y0_1_A1
    p_Y1_0_A0 <- 1 - p_Y1_1_A0
    
    p00 <- (1 - ps_eval) * p_Y0_0_A0;  p10 <- (1 - ps_eval) * p_Y0_1_A0
    p01 <- ps_eval       * p_Y0_0_A1;  p11 <- ps_eval       * p_Y0_1_A1
    
    ## alpha_0(y,x) — Gamma-free base odds ratio
    hat.OR.x.alpha0 <- safe_ratio(p11 * p00, p10 * p01)
    beta0           <- safe_ratio(p01, p00)
    
    ## Gamma-free conditional density ratio for DR term
    cdr_Y0 <- safe_ratio(ifelse(Y0 == 1, p_Y1_1_A0, p_Y1_0_A0),
                         ifelse(Y0 == 1, p_Y0_1_A0, p_Y0_0_A0))
    
  } else if (type == "poisson") {
    
    Y0r <- round(Y0); Y1r <- round(Y1); Ar <- round(A)
    df_all <- cbind(Y0 = Y0r, Y1 = Y1r, A = Ar, df_X)
    df_A0b <- df_all[Ar == 0, ]; df_A1b <- df_all[Ar == 1, ]
    
    fit_ps    <- stats::glm(A ~ ., data = df_all[, c("A", xnames)], family = stats::binomial())
    ps_eval   <- pmin(pmax(stats::predict(fit_ps, newdata = df_X, type = "response"), 1e-6), 1 - 1e-6)
    odds_X    <- safe_ratio(ps_eval, 1 - ps_eval)
    
    fit_Y0gA0 <- stats::glm(Y0 ~ ., data = df_A0b[, c("Y0", xnames)], family = stats::poisson())
    fit_Y0gA1 <- stats::glm(Y0 ~ ., data = df_A1b[, c("Y0", xnames)], family = stats::poisson())
    fit_Y1gA0 <- stats::glm(Y1 ~ ., data = df_A0b[, c("Y1", xnames)], family = stats::poisson())
    
    mu_Y0gA0 <- stats::predict(fit_Y0gA0, newdata = df_X, type = "response")
    mu_Y1gA0 <- stats::predict(fit_Y1gA0, newdata = df_X, type = "response")
    
    dr0 <- fit_density_ratio(Y0r[Ar == 1], X[Ar == 1, , drop = FALSE],
                             Y0r[Ar == 0], X[Ar == 0, , drop = FALSE])
    
    y_ref     <- round(stats::median(Y0r[Ar == 0]))
    r0_Y0     <- predict_density_ratio(dr0, Y0r,           X)
    r0_Y1     <- predict_density_ratio(dr0, Y1r,           X)
    r0_ref    <- predict_density_ratio(dr0, rep(y_ref, N),  X)
    
    ## alpha_0 (Gamma-free)
    alpha0_Y0 <- safe_ratio(r0_Y0, r0_ref)
    alpha0_Y1 <- safe_ratio(r0_Y1, r0_ref)
    beta0     <- r0_ref * odds_X
    
    ## Gamma-free grid summation
    y_max  <- max(stats::qpois(1 - 1e-6, lambda = max(mu_Y1gA0)), max(Y0r), max(Y1r))
    y_grid <- 0:y_max;  ng <- length(y_grid)
    
    X_rep_g    <- X[rep(seq_len(N), times = ng), , drop = FALSE]
    y_grid_rep <- rep(y_grid, each = N)
    alpha0_grid <- matrix(
      safe_ratio(predict_density_ratio(dr0, y_grid_rep,         X_rep_g),
                 predict_density_ratio(dr0, rep(y_ref, N * ng), X_rep_g)),
      N, ng)
    
    pmf_mat <- matrix(stats::dpois(rep(y_grid, each = N),
                            lambda = rep(mu_Y1gA0, times = ng)), N, ng)
    y_mat   <- matrix(y_grid, N, ng, byrow = TRUE)
    
    ## DR ratio (Gamma-free)
    dr_ratio <- safe_ratio(stats::dpois(Y0r, mu_Y1gA0), stats::dpois(Y0r, mu_Y0gA0))
    
  } else {
    stop("type must be 'continuous', 'binary', or 'poisson'")
  }
  
  ## ----------------------------------------------------------
  ## Step 1: Compute baseline mu(x) at Gamma = 1
  ## ----------------------------------------------------------
  
  .compute_eif_parametric <- function(type, Gamma, direction,
                                      mu_base_vec, mu_base_grid) {
    ## direction = "UB" or "LB" or "none" (Gamma=1)
    if (type == "continuous") {
      if (Gamma == 1) {
        alpha1_grid <- alpha0_grid
        alpha1_Y0   <- alpha0_Y0
        alpha1_Y1   <- alpha0_Y1
      } else {
        ## scale_grid: N x n_grid matrix; compare Y1_grid_mat against mu_base_grid
        if (direction == "UB") {
          scale_grid <- sens_scale_UB(Y1_grid_mat, mu_base_grid, Gamma)
          scale_Y0   <- sens_scale_UB(Y0, mu_base_vec, Gamma)
          scale_Y1   <- sens_scale_UB(Y1, mu_base_vec, Gamma)
        } else {
          scale_grid <- sens_scale_LB(Y1_grid_mat, mu_base_grid, Gamma)
          scale_Y0   <- sens_scale_LB(Y0, mu_base_vec, Gamma)
          scale_Y1   <- sens_scale_LB(Y1, mu_base_vec, Gamma)
        }
        alpha1_grid <- alpha0_grid * scale_grid
        alpha1_Y0   <- alpha0_Y0  * scale_Y0
        alpha1_Y1   <- alpha0_Y1  * scale_Y1
      }
      alpha1_grid[, !valid_grid] <- 0
      
      E_alpha1   <- rowSums(alpha1_grid * wt_mat)
      E_Y1alpha1 <- rowSums(Y1_grid_mat * alpha1_grid * wt_mat)
      
      hat.OR1           <- alpha1_Y1
      hat.BetaA1.plugin <- safe_ratio(odds_X, E_alpha1)
      hat.Mu1.plugin    <- safe_ratio(E_Y1alpha1, E_alpha1)
      hat.DR.plugin     <- hat.BetaA1.plugin *
        (A / pmax(beta0 * alpha0_Y0, 1e-6) + (1 - A)) * alpha1_Y0 * r1_Y0
      
    } else if (type == "binary") {
      ## Binary: y in {0,1}, mu_base_vec = P(Y1=1|A=0,X) under Gamma=1
      ## For y=1: compare 1 > mu => scale; for y=0: compare 0 > mu => scale
      if (Gamma == 1) {
        hat.OR.x.alpha1 <- hat.OR.x.alpha0
        alpha1_at_0     <- 1
      } else {
        ## alpha_0(1,x) = hat.OR.x.alpha0; alpha_0(0,x) = 1
        ## scale at y=1:
        if (direction == "UB") {
          s1 <- sens_scale_UB(rep(1, N), mu_base_vec, Gamma)
          s0 <- sens_scale_UB(rep(0, N), mu_base_vec, Gamma)
        } else {
          s1 <- sens_scale_LB(rep(1, N), mu_base_vec, Gamma)
          s0 <- sens_scale_LB(rep(0, N), mu_base_vec, Gamma)
        }
        hat.OR.x.alpha1 <- hat.OR.x.alpha0 * s1
        alpha1_at_0     <- 1 * s0  ## alpha_0(0,x)=1 scaled
      }
      
      E_alpha1   <- p_Y1_1_A0 * hat.OR.x.alpha1 + p_Y1_0_A0 * alpha1_at_0
      E_Y1alpha1 <- p_Y1_1_A0 * hat.OR.x.alpha1 * 1  ## Y1=1 contributes
      
      hat.OR1           <- ifelse(Y1 == 1, hat.OR.x.alpha1, alpha1_at_0)
      hat.BetaA1.plugin <- safe_ratio(odds_X, E_alpha1)
      hat.Mu1.plugin    <- safe_ratio(E_Y1alpha1, E_alpha1)
      hat.DR.plugin     <- hat.BetaA1.plugin *
        (A / pmax(beta0 * hat.OR.x.alpha0 ^ Y0, 1e-6) + (1 - A)) *
        ifelse(Y0 == 1, hat.OR.x.alpha1, alpha1_at_0) * cdr_Y0
      
    } else {  ## poisson
      if (Gamma == 1) {
        alpha1_grid <- alpha0_grid
        alpha1_Y0   <- alpha0_Y0
        alpha1_Y1   <- alpha0_Y1
      } else {
        if (direction == "UB") {
          scale_grid <- sens_scale_UB(y_mat, mu_base_grid, Gamma)
          scale_Y0   <- sens_scale_UB(Y0, mu_base_vec, Gamma)
          scale_Y1   <- sens_scale_UB(Y1, mu_base_vec, Gamma)
        } else {
          scale_grid <- sens_scale_LB(y_mat, mu_base_grid, Gamma)
          scale_Y0   <- sens_scale_LB(Y0, mu_base_vec, Gamma)
          scale_Y1   <- sens_scale_LB(Y1, mu_base_vec, Gamma)
        }
        alpha1_grid <- alpha0_grid * scale_grid
        alpha1_Y0   <- alpha0_Y0  * scale_Y0
        alpha1_Y1   <- alpha0_Y1  * scale_Y1
      }
      
      E_alpha1   <- rowSums(pmf_mat * alpha1_grid)
      E_Y1alpha1 <- rowSums(pmf_mat * alpha1_grid * y_mat)
      
      hat.OR1           <- alpha1_Y1
      hat.BetaA1.plugin <- safe_ratio(odds_X, E_alpha1)
      hat.Mu1.plugin    <- safe_ratio(E_Y1alpha1, E_alpha1)
      hat.DR.plugin     <- hat.BetaA1.plugin *
        (A / pmax(beta0 * alpha0_Y0, 1e-6) + (1 - A)) * alpha1_Y0 * dr_ratio
    }
    
    uncEIF <- A * Y1 / pr_A -
      EIF_Continuous(Y1, Y0, A, hat.OR1, hat.BetaA1.plugin,
                     hat.Mu1.plugin, hat.DR.plugin) / pr_A
    tau <- mean(uncEIF)
    se  <- stats::sd((uncEIF * pr_A - tau * A) / pr_A) / sqrt(N)
    list(ATT = tau, SE = se, Mu1 = hat.Mu1.plugin)
  }
  
  ## Baseline (Gamma = 1): get mu(x)
  base <- .compute_eif_parametric(type, Gamma = 1, direction = "none",
                                  mu_base_vec = NULL, mu_base_grid = NULL)
  mu_base_vec <- base$Mu1   # length N — per-unit baseline conditional mean
  
  ## Build mu_base in grid form for continuous / poisson
  if (type == "continuous") {
    mu_base_grid <- matrix(mu_base_vec, N, n_grid)
  } else if (type == "poisson") {
    mu_base_grid <- matrix(mu_base_vec, N, ng)
  } else {
    mu_base_grid <- NULL
  }
  
  ## ----------------------------------------------------------
  ## Sweep over Gamma_seq — compute UB and LB for each Gamma
  ## ----------------------------------------------------------
  
  results <- vector("list", length(Gamma_seq))
  
  for (gi in seq_along(Gamma_seq)) {
    Gamma <- Gamma_seq[gi]
    
    if (Gamma == 1) {
      results[[gi]] <- list(ATT_LB = base$ATT, SE_LB = base$SE,
                            ATT_UB = base$ATT, SE_UB = base$SE, Gamma = Gamma)
    } else {
      ## UB/LB for counterfactual mean => LB/UB for ATT
      res_UB <- .compute_eif_parametric(type, Gamma, "UB", mu_base_vec, mu_base_grid)
      res_LB <- .compute_eif_parametric(type, Gamma, "LB", mu_base_vec, mu_base_grid)
      results[[gi]] <- list(ATT_LB = res_UB$ATT, SE_LB = res_UB$SE,
                            ATT_UB = res_LB$ATT, SE_UB = res_LB$SE, Gamma = Gamma)
    }
  }
  
  ## Return
  effect_row <- data.frame(ATT = base$ATT, SE = base$SE)
  sens_df    <- data.frame(
    log_Gamma = log_Gamma_seq,
    ATT_LB    = sapply(results, `[[`, "ATT_LB"),
    SE_LB     = sapply(results, `[[`, "SE_LB"),
    ATT_UB    = sapply(results, `[[`, "ATT_UB"),
    SE_UB     = sapply(results, `[[`, "SE_UB")
  )
  return(list(Effect = effect_row, Sensitivity = sens_df))
}


## ============================================================
## NONPARAMETRIC SENSITIVITY ESTIMATOR (2-fold cross-fitting)
## ============================================================
