#' @keywords internal
UDID_Parametric_NoX <- function(Y0,
                                Y1,
                                A,
                                type          = "continuous",
                                log_Gamma_seq = 0) {
  
  Gamma_seq <- exp(log_Gamma_seq)
  N      <- length(Y0)
  idx0   <- which(A == 0)
  idx1   <- which(A == 1)
  pr_A   <- mean(A)
  
  ## ----------------------------------------------------------
  ## Fit all nuisance models ONCE, independent of Gamma
  ## ----------------------------------------------------------
  
  if (type == "continuous") {
    
    ## f1(Y1 | A=0) via Box-Cox + Gaussian (intercept-only)
    Y1_A0     <- Y1[idx0]
    shift_Y1  <- bc_find_shift(Y1_A0)
    Y1_A0_pos <- Y1_A0 + shift_Y1
    
    bc_Y1     <- MASS::boxcox(Y1_A0_pos ~ 1,
                              lambda = seq(-2, 2, by = 0.05), plotit = FALSE)
    lambda_Y1 <- bc_Y1$x[which.max(bc_Y1$y)]
    
    T_Y1_A0   <- bc_transform(Y1_A0_pos, lambda_Y1)
    mu1_bc    <- mean(T_Y1_A0)
    sigma2_f1 <- stats::var(T_Y1_A0)
    sig1_bc   <- sqrt(sigma2_f1)
    
    ## pr(A=1) — marginal, no X
    ps_eval <- rep(pr_A, N)
    odds_X  <- safe_ratio(ps_eval, 1 - ps_eval)
    
    ## Density ratio classifiers (Y only, Gamma-free)
    dr0 <- fit_density_ratio_NoX(Y0[idx1], Y0[idx0])
    dr1 <- fit_density_ratio_NoX(Y1[idx0], Y0[idx0])
    
    y_ref     <- stats::median(Y0[idx0])
    r0_Y0     <- predict_density_ratio_NoX(dr0, Y0)
    r0_Y1     <- predict_density_ratio_NoX(dr0, Y1)
    r0_ref    <- predict_density_ratio_NoX(dr0, rep(y_ref, N))
    
    ## alpha_0 (Gamma-free)
    alpha0_Y0 <- safe_ratio(r0_Y0, r0_ref)
    alpha0_Y1 <- safe_ratio(r0_Y1, r0_ref)
    beta0     <- r0_ref * odds_X
    r1_Y0     <- predict_density_ratio_NoX(dr1, Y0)
    
    ## Grid integration — equally-spaced grid in Box-Cox space
    n_grid    <- 501
    
    t_lo      <- mu1_bc - 4 * sig1_bc
    t_hi      <- mu1_bc + 4 * sig1_bc
    t_grid    <- seq(t_lo, t_hi, length.out = n_grid)
    dT        <- diff(t_grid)[1]
    
    Y1_grid_raw   <- bc_inverse(t_grid, lambda_Y1) - shift_Y1
    valid_grid    <- is.finite(Y1_grid_raw)
    Y1_grid_raw[!valid_grid] <- y_ref
    
    ## Density weights: phi((t_j - mu) / sig) / sig * dT  (same for all i)
    wt_vec <- stats::dnorm(t_grid, mean = mu1_bc, sd = sig1_bc) * dT
    wt_mat <- matrix(wt_vec, N, n_grid, byrow = TRUE)
    
    ## Evaluate alpha_0 on grid (no X dependence -> vectorize over grid only)
    alpha0_grid <- matrix(
      safe_ratio(rep(predict_density_ratio_NoX(dr0, Y1_grid_raw), each = N),
                 rep(predict_density_ratio_NoX(dr0, rep(y_ref, n_grid)), each = N)),
      N, n_grid)
    alpha0_grid[, !valid_grid] <- 0
    Y1_grid_mat <- matrix(Y1_grid_raw, N, n_grid, byrow = TRUE)
    Y1_grid_mat[, !valid_grid] <- 0
    
  } else if (type == "binary") {
    
    Y0r <- round(Y0); Y1r <- round(Y1); Ar <- round(A)
    
    ## Marginal probabilities (no X)
    ps_eval <- rep(pr_A, N)
    odds_X  <- safe_ratio(ps_eval, 1 - ps_eval)
    
    p_Y0_1_A0 <- rep(mean(Y0r[Ar == 0]), N)
    p_Y0_1_A1 <- rep(mean(Y0r[Ar == 1]), N)
    p_Y1_1_A0 <- rep(mean(Y1r[Ar == 0]), N)
    p_Y0_0_A0 <- 1 - p_Y0_1_A0;  p_Y0_0_A1 <- 1 - p_Y0_1_A1
    p_Y1_0_A0 <- 1 - p_Y1_1_A0
    
    p00 <- (1 - ps_eval) * p_Y0_0_A0;  p10 <- (1 - ps_eval) * p_Y0_1_A0
    p01 <- ps_eval       * p_Y0_0_A1;  p11 <- ps_eval       * p_Y0_1_A1
    
    hat.OR.x.alpha0 <- safe_ratio(p11 * p00, p10 * p01)
    beta0           <- safe_ratio(p01, p00)
    
    cdr_Y0 <- safe_ratio(ifelse(Y0 == 1, p_Y1_1_A0, p_Y1_0_A0),
                         ifelse(Y0 == 1, p_Y0_1_A0, p_Y0_0_A0))
    
  } else if (type == "poisson") {
    
    Y0r <- round(Y0); Y1r <- round(Y1); Ar <- round(A)
    
    ps_eval <- rep(pr_A, N)
    odds_X  <- safe_ratio(ps_eval, 1 - ps_eval)
    
    ## Marginal Poisson means (no X)
    mu_Y0gA0 <- rep(mean(Y0r[Ar == 0]), N)
    mu_Y1gA0 <- rep(mean(Y1r[Ar == 0]), N)
    
    dr0 <- fit_density_ratio_NoX(Y0r[Ar == 1], Y0r[Ar == 0])
    
    y_ref     <- round(stats::median(Y0r[Ar == 0]))
    r0_Y0     <- predict_density_ratio_NoX(dr0, Y0r)
    r0_Y1     <- predict_density_ratio_NoX(dr0, Y1r)
    r0_ref    <- predict_density_ratio_NoX(dr0, rep(y_ref, N))
    
    alpha0_Y0 <- safe_ratio(r0_Y0, r0_ref)
    alpha0_Y1 <- safe_ratio(r0_Y1, r0_ref)
    beta0     <- r0_ref * odds_X
    
    y_max  <- max(stats::qpois(1 - 1e-6, lambda = max(mu_Y1gA0)), max(Y0r), max(Y1r))
    y_grid <- 0:y_max;  ng <- length(y_grid)
    
    alpha0_grid <- matrix(
      safe_ratio(rep(predict_density_ratio_NoX(dr0, y_grid), each = N),
                 rep(predict_density_ratio_NoX(dr0, rep(y_ref, ng)), each = N)),
      N, ng)
    
    pmf_mat <- matrix(stats::dpois(rep(y_grid, each = N),
                            lambda = rep(mu_Y1gA0, times = ng)), N, ng)
    y_mat   <- matrix(y_grid, N, ng, byrow = TRUE)
    
    dr_ratio <- safe_ratio(stats::dpois(Y0r, mu_Y1gA0), stats::dpois(Y0r, mu_Y0gA0))
    
  } else {
    stop("type must be 'continuous', 'binary', or 'poisson'")
  }
  
  ## ----------------------------------------------------------
  ## Step 1: Compute baseline mu(x) at Gamma = 1
  ## ----------------------------------------------------------
  
  .compute_eif_noX <- function(type, Gamma, direction,
                               mu_base_vec, mu_base_grid) {
    if (type == "continuous") {
      if (Gamma == 1) {
        alpha1_grid <- alpha0_grid
        alpha1_Y0   <- alpha0_Y0
        alpha1_Y1   <- alpha0_Y1
      } else {
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
      if (Gamma == 1) {
        hat.OR.x.alpha1.loc <- hat.OR.x.alpha0
        alpha1_at_0         <- 1
      } else {
        if (direction == "UB") {
          s1 <- sens_scale_UB(rep(1, N), mu_base_vec, Gamma)
          s0 <- sens_scale_UB(rep(0, N), mu_base_vec, Gamma)
        } else {
          s1 <- sens_scale_LB(rep(1, N), mu_base_vec, Gamma)
          s0 <- sens_scale_LB(rep(0, N), mu_base_vec, Gamma)
        }
        hat.OR.x.alpha1.loc <- hat.OR.x.alpha0 * s1
        alpha1_at_0         <- 1 * s0
      }
      
      E_alpha1   <- p_Y1_1_A0 * hat.OR.x.alpha1.loc + p_Y1_0_A0 * alpha1_at_0
      E_Y1alpha1 <- p_Y1_1_A0 * hat.OR.x.alpha1.loc
      
      hat.OR1           <- ifelse(Y1 == 1, hat.OR.x.alpha1.loc, alpha1_at_0)
      hat.BetaA1.plugin <- safe_ratio(odds_X, E_alpha1)
      hat.Mu1.plugin    <- safe_ratio(E_Y1alpha1, E_alpha1)
      hat.DR.plugin     <- hat.BetaA1.plugin *
        (A / pmax(beta0 * hat.OR.x.alpha0 ^ Y0, 1e-6) + (1 - A)) *
        ifelse(Y0 == 1, hat.OR.x.alpha1.loc, alpha1_at_0) * cdr_Y0
      
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
  
  ## Baseline (Gamma = 1)
  base <- .compute_eif_noX(type, Gamma = 1, direction = "none",
                           mu_base_vec = NULL, mu_base_grid = NULL)
  mu_base_vec <- base$Mu1
  
  if (type == "continuous") {
    mu_base_grid <- matrix(mu_base_vec, N, n_grid)
  } else if (type == "poisson") {
    mu_base_grid <- matrix(mu_base_vec, N, ng)
  } else {
    mu_base_grid <- NULL
  }
  
  ## ----------------------------------------------------------
  ## Sweep over Gamma_seq — compute UB and LB
  ## ----------------------------------------------------------
  
  results <- vector("list", length(Gamma_seq))
  
  for (gi in seq_along(Gamma_seq)) {
    Gamma <- Gamma_seq[gi]
    
    if (Gamma == 1) {
      results[[gi]] <- list(ATT_LB = base$ATT, SE_LB = base$SE,
                            ATT_UB = base$ATT, SE_UB = base$SE, Gamma = Gamma)
    } else {
      ## UB/LB for counterfactual mean => LB/UB for ATT
      res_UB <- .compute_eif_noX(type, Gamma, "UB", mu_base_vec, mu_base_grid)
      res_LB <- .compute_eif_noX(type, Gamma, "LB", mu_base_vec, mu_base_grid)
      results[[gi]] <- list(ATT_LB = res_UB$ATT, SE_LB = res_UB$SE,
                            ATT_UB = res_LB$ATT, SE_UB = res_LB$SE, Gamma = Gamma)
    }
  }
  
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
## NONPARAMETRIC ESTIMATOR — NO COVARIATES (2-fold cross-fitting)
## ============================================================

