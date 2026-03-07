#' @keywords internal
UDID_Nonparametric_NoX <- function(Y0,
                                   Y1,
                                   A,
                                   type                  = "continuous",
                                   log_Gamma_seq         = 0,
                                   seed                  = 42,
                                   SL.list               = c(1),
                                   hyperparameter        = "fast",
                                   density.report        = FALSE,
                                   hyperparameter.report = FALSE,
                                   verbose               = FALSE) {
  Gamma_seq <- exp(log_Gamma_seq)
  N  <- length(Y0)
  nG <- length(Gamma_seq)
  
  ## Stratified 2-fold split (same as UDID_Nonparametric)
  set.seed(seed)
  idx1     <- sample(which(A == 1), floor(sum(A == 1) / 2))
  idx0     <- sample(which(A == 0), ceiling(sum(A == 0) / 2))
  SS.Index <- list(c(idx1, idx0), setdiff(seq_len(N), c(idx1, idx0)))
  
  uncEIF_UB_full   <- replicate(nG, numeric(N), simplify = FALSE)
  uncEIF_LB_full   <- replicate(nG, numeric(N), simplify = FALSE)
  uncEIF_base_full <- numeric(N)   ## EIF at Gamma = 1, always accumulated
  
  ## Density accumulators — continuous: N x 501; binary: N x 2
  if (density.report) {
    Num.Grid <- if (type == "continuous") 501L else 2L
    density_full <- list(
      fY0A0  = matrix(NA_real_, N, Num.Grid),
      fY0A1  = matrix(NA_real_, N, Num.Grid),
      fY1A0  = matrix(NA_real_, N, Num.Grid),
      fY1A1  = matrix(NA_real_, N, Num.Grid),
      Y.Grid = NULL
    )
  }
  
  ## Hyperparameter report: (continuous only) bandwidth info
  if (hyperparameter.report && type == "continuous") {
    fold_bw_info <- vector("list", 2)
  }
  
  for (ss in 1:2) {
    train <- SS.Index[[ss]]; eval <- SS.Index[[3 - ss]]
    fold_results <- .UDID_fold_NoX(
      Y0[train], Y1[train], A[train],
      Y0[eval],  Y1[eval],  A[eval],
      type, SL.list, Gamma_seq, seed = seed * 10L + ss, fold_id = ss,
      hyperparameter        = hyperparameter,
      density.report        = density.report,
      hyperparameter.report = hyperparameter.report,
      verbose               = verbose)
    for (gi in seq_len(nG)) {
      ## UB/LB for counterfactual mean => LB/UB for ATT
      uncEIF_UB_full[[gi]][eval] <- fold_results[[gi]]$LB
      uncEIF_LB_full[[gi]][eval] <- fold_results[[gi]]$UB
    }
    uncEIF_base_full[eval] <- fold_results$base_eif
    if (density.report && !is.null(fold_results$fY0A0)) {
      density_full$fY0A0[eval, ] <- fold_results$fY0A0
      density_full$fY0A1[eval, ] <- fold_results$fY0A1
      density_full$fY1A0[eval, ] <- fold_results$fY1A0
      density_full$fY1A1[eval, ] <- fold_results$fY1A1
      density_full$Y.Grid        <- fold_results$Y.Grid
    }
    if (hyperparameter.report && type == "continuous") {
      fold_bw_info[[ss]] <- fold_results$bw_info
    }
  }
  
  pr_A    <- mean(A)
  results <- vector("list", nG)
  
  base_gi <- which(Gamma_seq == 1)
  
  for (gi in seq_len(nG)) {
    ## UB
    uncEIF_UB  <- uncEIF_UB_full[[gi]] / pr_A
    tau_UB     <- mean(uncEIF_UB)
    cEIF_UB    <- uncEIF_UB - A * tau_UB / pr_A
    SE_UB      <- stats::sd(cEIF_UB) / sqrt(N)
    set.seed(seed + gi)
    SE.Boot_UB <- stats::sd(Mboot(uncEIF_UB, N, NumBoot = 10000))
    
    ## LB
    uncEIF_LB  <- uncEIF_LB_full[[gi]] / pr_A
    tau_LB     <- mean(uncEIF_LB)
    cEIF_LB    <- uncEIF_LB - A * tau_LB / pr_A
    SE_LB      <- stats::sd(cEIF_LB) / sqrt(N)
    set.seed(seed + gi + 1000L)
    SE.Boot_LB <- stats::sd(Mboot(uncEIF_LB, N, NumBoot = 10000))
    
    results[[gi]] <- list(ATT_LB = tau_LB, SE_LB = SE_LB, SE.Boot_LB = SE.Boot_LB,
                          ATT_UB = tau_UB, SE_UB = SE_UB, SE.Boot_UB = SE.Boot_UB,
                          Gamma  = Gamma_seq[gi])
  }
  
  if (length(base_gi) > 0) {
    base_att <- results[[base_gi[1]]]$ATT_UB
    base_se  <- results[[base_gi[1]]]$SE_UB
    base_bse <- results[[base_gi[1]]]$SE.Boot_UB
  } else {
    base_att <- results[[1]]$ATT_UB
    base_se  <- results[[1]]$SE_UB
    base_bse <- results[[1]]$SE.Boot_UB
  }
  effect_row <- data.frame(ATT = base_att, SE = base_se, Boot.SE = base_bse)
  
  sens_df <- data.frame(
    log_Gamma  = log_Gamma_seq,
    ATT_LB     = sapply(results, `[[`, "ATT_LB"),
    SE_LB      = sapply(results, `[[`, "SE_LB"),
    SE.Boot_LB = sapply(results, `[[`, "SE.Boot_LB"),
    ATT_UB     = sapply(results, `[[`, "ATT_UB"),
    SE_UB      = sapply(results, `[[`, "SE_UB"),
    SE.Boot_UB = sapply(results, `[[`, "SE.Boot_UB")
  )
  
  ## EIF at Gamma = 1 (always reported)
  eif_unc <- uncEIF_base_full / pr_A
  eif_tau <- mean(eif_unc)
  eif_cen <- eif_unc - A * eif_tau / pr_A
  eif_df  <- data.frame(
    `Uncentered EIF` = eif_unc,
    `Centered EIF`   = eif_cen,
    check.names = FALSE
  )
  
  RESULT <- list(Effect = effect_row, EIF = eif_df, Sensitivity = sens_df)
  
  ## Hyperparameter report: SL_Library = NULL (no SuperLearner in NoX),
  ## selected_bws for continuous case
  if (hyperparameter.report) {
    RESULT$SL_Library <- NULL
    
    if (type == "continuous") {
      bw_i1 <- fold_bw_info[[1]]
      bw_i2 <- fold_bw_info[[2]]
      bw_table <- data.frame(
        `Fold 1: f(Y0|A=0)` = bw_i1$fY0A0$ybw,
        `Fold 2: f(Y0|A=0)` = bw_i2$fY0A0$ybw,
        `Fold 1: f(Y0|A=1)/f(Y0|A=0)` = bw_i1$DR_Y0A1_Y0A0$sigma,
        `Fold 2: f(Y0|A=1)/f(Y0|A=0)` = bw_i2$DR_Y0A1_Y0A0$sigma,
        `Fold 1: f(Y1|A=0)/f(Y0|A=0)` = bw_i1$DR_Y1A0_Y0A0$sigma,
        `Fold 2: f(Y1|A=0)/f(Y0|A=0)` = bw_i2$DR_Y1A0_Y0A0$sigma,
        row.names   = "Y",
        check.names = FALSE
      )
      RESULT$selected_bws <- bw_table
    }
  }
  
  if (density.report) {
    RESULT$Density <- density_full
  }
  return(RESULT)
}


## ---- Fold worker (no covariates) ----------------------------

.UDID_fold_NoX <- function(Y0, Y1, A,
                           Y0.Eval, Y1.Eval, A.Eval,
                           type, SL.list, Gamma_seq, seed, fold_id = NULL,
                           hyperparameter        = "fast",
                           density.report        = FALSE,
                           hyperparameter.report = FALSE,
                           verbose               = FALSE) {
  
  N.Test           <- length(Y0.Eval)
  Num.Y.Grid.Basis <- 501
  
  Y.Grid.Basis <- seq(
    min(c(Y0, Y1[A == 0])) - max(stats::sd(Y0), stats::sd(Y1[A == 0])),
    max(c(Y0, Y1[A == 0])) + max(stats::sd(Y0), stats::sd(Y1[A == 0])),
    length = Num.Y.Grid.Basis
  )
  dY        <- diff(Y.Grid.Basis)[1]
  grid_cols <- 4:(3 + Num.Y.Grid.Basis)
  
  ## ----------------------------------------------------------
  ## Fit all nuisance models ONCE, independent of Gamma
  ## ----------------------------------------------------------
  
  if (type == "continuous") {
    
    ## Evaluation y-vectors: (Y0.Eval, Y1.Eval, y_ref, Y.Grid.Basis)
    eval_y_vec <- c(Y0.Eval, Y1.Eval,
                    rep(stats::median(Y0[A == 0]), N.Test),
                    rep(Y.Grid.Basis, each = N.Test))
    
    ## KDE of Y0 | A=0 (marginal, no X conditioning)
    set.seed(seed)
    if (verbose) cat(sprintf("Fold %d: estimating the density f(Y_{t=0}|A=0)\n", fold_id))
    if(hyperparameter=="fast"){
      npbw.Y0A0    <- np::npcdensbw(xdat = factor(rep(1, sum(A == 0))), ydat = Y0[A == 0], bwmethod="normal-reference")
    } else {
      npbw.Y0A0    <- np::npcdensbw(xdat = factor(rep(1, sum(A == 0))), ydat = Y0[A == 0])
    }
    
    kde_vals     <- stats::fitted(np::npcdens(npbw.Y0A0,
                                              exdat = factor(rep(1, length(eval_y_vec))),
                                              eydat = eval_y_vec))
    kde_mat      <- matrix(kde_vals, N.Test)
    Raw.Y0A0.KDE <- kde_mat / (rowSums(kde_mat[, grid_cols]) * dY)
    
    ## Density ratio Y0|A=1 vs Y0|A=0 (Y only, no X)
    
    if (hyperparameter=="fast"){
      DR.sigma <- "median"
    } else {
      DR.sigma <- "auto"
    }
    
    set.seed(seed + 1)
    if (verbose) cat(sprintf("Fold %d: estimating the density ratio f(Y_{t=0}|A=1)/f(Y_{t=0}|A=0)\n", fold_id))
    ## Standardise Y for KLIEP; use Y0 mean/sd for all Y variables
    mu_Y <- mean(Y0);  sd_Y <- stats::sd(Y0)
    Y0s  <- (Y0 - mu_Y) / sd_Y
    Y1s  <- (Y1 - mu_Y) / sd_Y
    eval_y_vec_s <- (eval_y_vec - mu_Y) / sd_Y
    
    DR.Y0A1.Y0A0 <- KLIEP_bw(matrix(Y0s[A == 1]),
                             matrix(Y0s[A == 0]),
                             n_y_cols = 1,
                             sigma=DR.sigma,
                             verbose = verbose)
    Raw.Y0A1.DR  <- Raw.Y0A0.KDE *
      matrix(compute_density_ratio_Kernel(matrix(eval_y_vec_s), DR.Y0A1.Y0A0), N.Test)
    Raw.Y0A1.DR  <- Raw.Y0A1.DR / (rowSums(Raw.Y0A1.DR[, grid_cols]) * dY)
    
    ## Density ratio Y1|A=0 vs Y0|A=0 (Y only, no X)
    set.seed(seed + 2)
    if (verbose) cat(sprintf("Fold %d: estimating the density ratio f(Y_{t=1}|A=0)/f(Y_{t=0}|A=0)\n", fold_id))
    DR.Y1A0.Y0A0 <- KLIEP_bw(matrix(Y1s[A == 0]),
                             matrix(Y0s[A == 0]),
                             n_y_cols=1,
                             sigma=DR.sigma,
                             verbose = verbose)
    Raw.Y1A0.DR  <- Raw.Y0A0.KDE *
      matrix(compute_density_ratio_Kernel(matrix(eval_y_vec_s), DR.Y1A0.Y0A0), N.Test)
    Raw.Y1A0.DR  <- Raw.Y1A0.DR / (rowSums(Raw.Y1A0.DR[, grid_cols]) * dY)
    
    col_Y0 <- 1L; col_Y1 <- 2L; col_base <- 3L
    
    ## P(A=1) — marginal, no X
    if (verbose) cat(sprintf("Fold %d: estimating the propensity score Pr(A)\n", fold_id))
    Est.Marginal.ProbA <- rep(mean(A), N.Test)
    
    ## Capture KDE and KLIEP bandwidths (continuous only)
    if (hyperparameter.report) {
      bw_info <- list(
        fY0A0        = list(ybw = as.numeric(npbw.Y0A0$ybw)),
        DR_Y0A1_Y0A0 = list(sigma = DR.Y0A1.Y0A0$kernel_info$sigma),
        DR_Y1A0_Y0A0 = list(sigma = DR.Y1A0.Y0A0$kernel_info$sigma)
      )
    }
    
    if (verbose) cat(sprintf("Fold %d: calculating the efficient influence function\n", fold_id))
    Joint.Prob.Y0A.EvalY0 <- cbind(
      Raw.Y0A1.DR[, col_Y0]    * Est.Marginal.ProbA,
      Raw.Y0A1.DR[, col_base]  * Est.Marginal.ProbA,
      Raw.Y0A0.KDE[, col_Y0]   * (1 - Est.Marginal.ProbA),
      Raw.Y0A0.KDE[, col_base] * (1 - Est.Marginal.ProbA)
    )
    Joint.Prob.Y0A.EvalY1 <- cbind(
      Raw.Y0A1.DR[, col_Y1]    * Est.Marginal.ProbA,
      Raw.Y0A1.DR[, col_base]  * Est.Marginal.ProbA,
      Raw.Y0A0.KDE[, col_Y1]   * (1 - Est.Marginal.ProbA),
      Raw.Y0A0.KDE[, col_base] * (1 - Est.Marginal.ProbA)
    )
    
    hat.OR0.alpha0 <- Truncate.Function(
      Joint.Prob.Y0A.EvalY0[, 1] * Joint.Prob.Y0A.EvalY0[, 4] /
        (Joint.Prob.Y0A.EvalY0[, 2] * Joint.Prob.Y0A.EvalY0[, 3]), 0, 1e10)
    hat.OR1.alpha0 <- Truncate.Function(
      Joint.Prob.Y0A.EvalY1[, 1] * Joint.Prob.Y0A.EvalY1[, 4] /
        (Joint.Prob.Y0A.EvalY1[, 2] * Joint.Prob.Y0A.EvalY1[, 3]), 0, 1e10)
    
    hat.BetaA0 <- Joint.Prob.Y0A.EvalY0[, 2] / Joint.Prob.Y0A.EvalY0[, 4]
    DR.base    <- (Raw.Y1A0.DR[, col_Y0] / Raw.Y0A0.KDE[, col_Y0])
    
    OR.Grid.alpha0 <- Truncate.Function(
      Raw.Y0A1.DR[, grid_cols] *
        matrix(Raw.Y0A0.KDE[, col_base], N.Test, Num.Y.Grid.Basis) /
        (matrix(Raw.Y0A1.DR[, col_base], N.Test, Num.Y.Grid.Basis) *
           Raw.Y0A0.KDE[, grid_cols]),
      0, 1e10)
    Cond.Density <- Raw.Y1A0.DR[, grid_cols] * dY
    
  } else {
    ## Binary outcome — no covariates
    Data.Reset    <- data.frame(Y0 = round(Y0), Y1 = round(Y1), A = round(A))
    Data.Reset.A0 <- Data.Reset[Data.Reset$A == 0, ]
    Data.Reset.A1 <- Data.Reset[Data.Reset$A == 1, ]
    
    ## Marginal probabilities (no X, no MySL needed)
    if (verbose) cat(sprintf("Fold %d: estimating the propensity score Pr(A)\n", fold_id))
    Est.PrA1     <- rep(mean(A), N.Test)
    if (verbose) cat(sprintf("Fold %d: estimating the outcome regression Pr(Y_{t=0}|A=1)\n", fold_id))
    Est.PrY01gA1 <- rep(mean(round(Y0[A == 1])), N.Test)
    if (verbose) cat(sprintf("Fold %d: estimating the outcome regression Pr(Y_{t=0}|A=0)\n", fold_id))
    Est.PrY01gA0 <- rep(mean(round(Y0[A == 0])), N.Test)
    if (verbose) cat(sprintf("Fold %d: estimating the outcome regression Pr(Y_{t=1}|A=0)\n", fold_id))
    Est.PrY11gA0 <- rep(mean(round(Y1[A == 0])), N.Test)
    
    if (verbose) cat(sprintf("Fold %d: calculating the efficient influence function\n", fold_id))
    Prob.Predict <- cbind(
      (1 - Est.PrA1) * (1 - Est.PrY01gA0),
      (1 - Est.PrA1) * Est.PrY01gA0,
      Est.PrA1       * (1 - Est.PrY01gA1),
      Est.PrA1       * Est.PrY01gA1,
      (1 - Est.PrA1) * (1 - Est.PrY11gA0),
      (1 - Est.PrA1) * Est.PrY11gA0
    )
    
    hat.OR.x.alpha0    <- Prob.Predict[,4]*Prob.Predict[,1] / (Prob.Predict[,2]*Prob.Predict[,3])
    Est.Marginal.ProbA <- Prob.Predict[,3] + Prob.Predict[,4]
    Cond.Density       <- Prob.Predict[,6] / (Prob.Predict[,5] + Prob.Predict[,6])
    
    DR.base    <- (Prob.Predict[,5]*(1-Y0.Eval) + Prob.Predict[,6]*Y0.Eval) /
      (Prob.Predict[,1]*(1-Y0.Eval) + Prob.Predict[,2]*Y0.Eval)
    beta0.eval <- Prob.Predict[,3] / Prob.Predict[,1]
    
    ## Binary conditional densities: N x 2 matrices (col 1 = Y=0, col 2 = Y=1)
    ## fY1A1 is derived via the UDID OR assumption at Gamma = 1
    if (density.report) {
      bin_fY0A0      <- cbind(1 - Est.PrY01gA0, Est.PrY01gA0)
      bin_fY0A1      <- cbind(1 - Est.PrY01gA1, Est.PrY01gA1)
      bin_fY1A0      <- cbind(1 - Est.PrY11gA0, Est.PrY11gA0)
      denom_Y1A1     <- (1 - Cond.Density) + hat.OR.x.alpha0 * Cond.Density
      bin_fY1A1      <- cbind((1 - Cond.Density)              / denom_Y1A1,
                              hat.OR.x.alpha0 * Cond.Density / denom_Y1A1)
    }
  }
  
  ## ----------------------------------------------------------
  ## Sweep over Gamma_seq — compute UB and LB
  ## ----------------------------------------------------------
  
  pr_A.eval <- mean(A.Eval)
  if (type == "continuous") {
    Y.Grid.Mat <- matrix(Y.Grid.Basis, N.Test, Num.Y.Grid.Basis, byrow = TRUE)
  }
  
  .fold_eif_noX <- function(Gamma, direction) {
    if (type == "continuous") {
      if (Gamma == 1) {
        OR.Grid.alpha1 <- OR.Grid.alpha0
        hat.OR0.alpha1 <- hat.OR0.alpha0
        hat.OR1.alpha1 <- hat.OR1.alpha0
      } else {
        if (direction == "UB") {
          scale_grid <- sens_scale_UB(Y.Grid.Mat, mu_base_grid, Gamma)
          scale_Y0   <- sens_scale_UB(Y0.Eval, mu_base_vec, Gamma)
          scale_Y1   <- sens_scale_UB(Y1.Eval, mu_base_vec, Gamma)
        } else {
          scale_grid <- sens_scale_LB(Y.Grid.Mat, mu_base_grid, Gamma)
          scale_Y0   <- sens_scale_LB(Y0.Eval, mu_base_vec, Gamma)
          scale_Y1   <- sens_scale_LB(Y1.Eval, mu_base_vec, Gamma)
        }
        OR.Grid.alpha1 <- OR.Grid.alpha0 * scale_grid
        hat.OR0.alpha1 <- hat.OR0.alpha0 * scale_Y0
        hat.OR1.alpha1 <- hat.OR1.alpha0 * scale_Y1
      }
      
      E_alpha1   <- rowSums(OR.Grid.alpha1 * Cond.Density)
      E_Y1alpha1 <- rowSums(Y.Grid.Mat * OR.Grid.alpha1 * Cond.Density)
      
      hat.OR1           <- hat.OR1.alpha1
      hat.Mu1.plugin    <- E_Y1alpha1 / E_alpha1
      hat.BetaA1.plugin <- Est.Marginal.ProbA / (1 - Est.Marginal.ProbA) / E_alpha1
      hat.DR.withoutbA1 <- (A.Eval / (hat.BetaA0 * hat.OR0.alpha0) + (1 - A.Eval)) * hat.OR0.alpha1 * DR.base
      hat.DR.plugin     <- hat.DR.withoutbA1 * hat.BetaA1.plugin
      
    } else {
      ## Binary
      if (Gamma == 1) {
        hat.OR.x.alpha1.loc <- hat.OR.x.alpha0
        alpha1_at_0         <- 1
      } else {
        if (direction == "UB") {
          s1 <- sens_scale_UB(rep(1, N.Test), mu_base_vec, Gamma)
          s0 <- sens_scale_UB(rep(0, N.Test), mu_base_vec, Gamma)
        } else {
          s1 <- sens_scale_LB(rep(1, N.Test), mu_base_vec, Gamma)
          s0 <- sens_scale_LB(rep(0, N.Test), mu_base_vec, Gamma)
        }
        hat.OR.x.alpha1.loc <- hat.OR.x.alpha0 * s1
        alpha1_at_0         <- 1 * s0
      }
      
      E_alpha1 <- Cond.Density * hat.OR.x.alpha1.loc + (1 - Cond.Density) * alpha1_at_0
      
      hat.OR1           <- ifelse(Y1.Eval == 1, hat.OR.x.alpha1.loc, alpha1_at_0)
      hat.Mu1.plugin    <- hat.OR.x.alpha1.loc * Cond.Density / E_alpha1
      hat.BetaA1.plugin <- Est.Marginal.ProbA / (1 - Est.Marginal.ProbA) / E_alpha1
      hat.DR.plugin     <- (A.Eval / (beta0.eval * hat.OR.x.alpha0 ^ Y0.Eval) + (1 - A.Eval)) *
        ifelse(Y0.Eval == 1, hat.OR.x.alpha1.loc, alpha1_at_0) * DR.base * hat.BetaA1.plugin
    }
    
    uncEIF <- A.Eval * Y1.Eval / pr_A.eval -
      EIF_Continuous(Y1.Eval, Y0.Eval, A.Eval,
                     hat.OR1, hat.BetaA1.plugin, hat.Mu1.plugin,
                     hat.DR.plugin) / pr_A.eval
    list(eif = uncEIF * pr_A.eval, Mu1 = hat.Mu1.plugin)
  }
  
  ## Baseline at Gamma = 1
  base_res    <- .fold_eif_noX(Gamma = 1, direction = "none")
  mu_base_vec <- base_res$Mu1
  if (type == "continuous") {
    mu_base_grid <- matrix(mu_base_vec, N.Test, Num.Y.Grid.Basis)
  } else {
    mu_base_grid <- NULL
  }
  
  results <- vector("list", length(Gamma_seq))
  
  for (gi in seq_along(Gamma_seq)) {
    Gamma <- Gamma_seq[gi]
    
    if (Gamma == 1) {
      results[[gi]] <- list(UB = base_res$eif, LB = base_res$eif)
    } else {
      res_UB <- .fold_eif_noX(Gamma, "UB")
      res_LB <- .fold_eif_noX(Gamma, "LB")
      results[[gi]] <- list(UB = res_UB$eif, LB = res_LB$eif)
    }
  }
  
  ## Always store EIF at Gamma = 1 for the caller
  results$base_eif <- base_res$eif
  
  ## Density output
  if (density.report && type == "continuous") {
    results$fY0A0      <- Raw.Y0A0.KDE[, grid_cols]
    results$fY0A1      <- Raw.Y0A1.DR[, grid_cols]
    results$fY1A0      <- Raw.Y1A0.DR[, grid_cols]
    fY1A1              <- OR.Grid.alpha0 * Raw.Y1A0.DR[, grid_cols] * dY
    results$fY1A1      <- fY1A1 / matrix(rowSums(fY1A1), N.Test, Num.Y.Grid.Basis)
    results$Y.Grid     <- Y.Grid.Basis
  } else if (density.report && type == "binary") {
    results$fY0A0  <- bin_fY0A0
    results$fY0A1  <- bin_fY0A1
    results$fY1A0  <- bin_fY1A0
    results$fY1A1  <- bin_fY1A1
    results$Y.Grid <- c(0L, 1L)
  }
  
  ## Hyperparameter report output
  if (hyperparameter.report && type == "continuous") {
    results$bw_info <- bw_info
  }
  
  results
}
