#' UDID Parametric Sensitivity Parameter
#'
#' Estimate the range of the sensitivity parameter from observed data
#' using parametric density ratio classifiers.
#'
#' @param Yn1 Numeric vector of negative-one-period outcomes.
#' @param Y0 Numeric vector of pre-treatment outcomes.
#' @param A Binary treatment indicator (0/1).
#' @param X Optional covariate matrix.
#' @param type Outcome type: \code{"continuous"}, \code{"binary"}, or \code{"poisson"}.
#'
#' @return A named numeric vector of quantiles of the sensitivity parameter ratio.
#' @export
UDID_Parametric_SensPara <- function(Yn1,
                                     Y0,
                                     A,
                                     X    = NULL,
                                     type = "continuous",
                                     quantile.range = seq(0,1,by=0.005)) {
  
  ## Dispatch to no-covariate version if X is NULL
  if (is.null(X)) {
    return(UDID_Parametric_NoX_SensPara(Yn1 = Yn1, Y0 = Y0, A = A,
                                        type = type, quantile.range = quantile.range))
  }
  
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
  
  ## Density ratio classifiers (Gamma-free)
  dr0 <- fit_density_ratio(Y0[idx1], X[idx1, , drop = FALSE],
                           Y0[idx0], X[idx0, , drop = FALSE])
  drn1 <- fit_density_ratio(Yn1[idx1], X[idx1, , drop = FALSE],
                            Yn1[idx0], X[idx0, , drop = FALSE])
  
  if (type == "continuous") {
    
    y_ref     <- stats::median(Y0[idx0])
    
    ## ----- Grid setup to match nonparametric -----
    eval_y_vec <- c(Y0, Yn1, rep(y_ref, N))
    eval_x_mat <- X[rep(seq_len(N), times = 3), , drop = FALSE] 
    
    r0  <- matrix(predict_density_ratio(dr0,  eval_y_vec, eval_x_mat), N)
    rn1 <- matrix(predict_density_ratio(drn1, eval_y_vec, eval_x_mat), N)
    
    ## OR_Y0(y,x) = f(y|A=1,X)*f(y_ref|A=0,X) / (f(y_ref|A=1,X)*f(y|A=0,X))
    OR.Y0.alpha0  <- r0[,1] / r0[,3]
    OR.Yn1.alpha0 <- r0[,2] / r0[,3]
    
    OR.Y0.alphan1  <- rn1[,1] / rn1[,3]
    OR.Yn1.alphan1 <- rn1[,2] / rn1[,3]
    
    range_alpha <- c( OR.Y0.alpha0[idx1]  / OR.Y0.alphan1[idx1],
                      OR.Yn1.alpha0[idx1] / OR.Yn1.alphan1[idx1] )
    
    return( stats::quantile(range_alpha, quantile.range) )
    
  } else if (type == "binary"){
    
    y_ref     <- 0
    y_grid    <- c(0,1)
    ng <- length(y_grid)
    
    X_rep_g    <- X[rep(seq_len(N), times = ng), , drop = FALSE]
    y_grid_rep <- rep(y_grid, each = N)
    
    r0        <- predict_density_ratio(dr0,  y_grid_rep, X_rep_g)
    rn1       <- predict_density_ratio(drn1, y_grid_rep, X_rep_g)
    r0_ref    <- predict_density_ratio(dr0,  rep(y_ref, N * ng), X_rep_g)
    rn1_ref   <- predict_density_ratio(drn1, rep(y_ref, N * ng), X_rep_g)
    
    ## alpha_0 (Gamma-free)
    alpha0_grid  <- matrix(safe_ratio(r0, r0_ref),N,ng)
    alphan1_grid <- matrix(safe_ratio(rn1, rn1_ref),N,ng)
    range_alpha <- as.numeric(alpha0_grid/alphan1_grid)
    
    return( stats::quantile(range_alpha,quantile.range) )
    
    
  } else if (type == "poisson"){
    
    y_ref     <- stats::median(Y0[idx0])
    y_grid    <- c(0,max(c(Y0[idx1],Yn1[idx1])))
    ng <- length(y_grid)
    
    X_rep_g    <- X[rep(seq_len(N), times = ng), , drop = FALSE]
    y_grid_rep <- rep(y_grid, each = N)
    
    r0        <- predict_density_ratio(dr0,  y_grid_rep, X_rep_g)
    rn1       <- predict_density_ratio(drn1, y_grid_rep, X_rep_g)
    r0_ref    <- predict_density_ratio(dr0,  rep(y_ref, N * ng), X_rep_g)
    rn1_ref   <- predict_density_ratio(drn1, rep(y_ref, N * ng), X_rep_g)
    
    ## alpha_0 (Gamma-free)
    alpha0_grid  <- matrix(safe_ratio(r0, r0_ref),N,ng)
    alphan1_grid <- matrix(safe_ratio(rn1, rn1_ref),N,ng)
    range_alpha  <- as.numeric(alpha0_grid / alphan1_grid)
    
    return( stats::quantile(range_alpha,quantile.range) )
    
    
  } else {
    stop("type must be 'continuous', 'binary', or 'poisson'")
  }
  
}






## ============================================================
## NONPARAMETRIC — SENSITIVITY PARAMETER
## ============================================================

#' UDID Nonparametric Sensitivity Parameter
#'
#' Estimate the range of the sensitivity parameter from observed data
#' using nonparametric Super Learner models.
#'
#' @param Yn1 Numeric vector of negative-one-period outcomes.
#' @param Y0 Numeric vector of pre-treatment outcomes.
#' @param A Binary treatment indicator (0/1).
#' @param X Optional covariate matrix.
#' @param type Outcome type: \code{"continuous"} or \code{"binary"}.
#' @param seed Random seed.
#' @param hyperparameter Hyperparameter tuning \code{"fast"} or \code{"slow"}.
#' @param SL.hpara Super Learner hyperparameters.
#'
#' @return A named numeric vector of quantiles of the sensitivity parameter ratio.
#' @export
UDID_Nonparametric_SensPara <- function(Yn1,
                                        Y0,
                                        A,
                                        X         = NULL,
                                        type      = "continuous",
                                        quantile.range = seq(0,1,by=0.005),
                                        seed      = 42,
                                        hyperparameter = "fast",
                                        SL.hpara  = list(SLL      = c(1),
                                                         MLPL     = c(2, 4),
                                                         MTRY     = c(1, 2, 3),
                                                         NMN      = 25,
                                                         MLPdecay = c(0.1, 0.001))) {
  
  ## Dispatch to no-covariate version if X is NULL
  if (is.null(X)) {
    return(UDID_Nonparametric_NoX_SensPara(Yn1 = Yn1, Y0 = Y0, A = A,
                                           type = type, seed = seed,
                                           hyperparameter = hyperparameter,
                                           quantile.range = quantile.range,
                                           SL.hpara = SL.hpara))
  }
  
  N      <- length(Y0)
  X      <- as.matrix(X)
  d      <- ncol(X)
  xnames <- paste0("X", sprintf("%05d", seq_len(d)))
  df_X   <- stats::setNames(as.data.frame(X), xnames)
  idx0   <- which(A == 0)
  idx1   <- which(A == 1)
  pr_A   <- mean(A)
  
  if (type == "continuous") {
    
    ## ----- Density ratios via KLIEP_bw (per-dimension bandwidth) -----
    ## Standardise Y and X; use Y0 mean/sd for all Y variables
    mu_Y   <- mean(Y0);     sd_Y   <- stats::sd(Y0)
    mu_X   <- colMeans(X);  sd_X <- apply(X, 2, sd)
    sd_X[sd_X == 0] <- 1                          # guard constant columns
    
    Y0s  <- (Y0  - mu_Y)  / sd_Y
    Yn1s <- (Yn1 - mu_Y)  / sd_Y
    Xs   <- sweep(sweep(X, 2, mu_X), 2, sd_X, "/")
    
    y_ref     <- stats::median(Y0[idx0])
    
    ## ----- Grid setup -----
    eval_y_vec <- c(Y0, Yn1, rep(y_ref, N))
    eval_x_mat <- X[rep(seq_len(N), times = 3), , drop = FALSE] 
    
    eval_y_vec_s <- (eval_y_vec - mu_Y) / sd_Y
    eval_x_mat_s <- sweep(sweep(eval_x_mat, 2, mu_X), 2, sd_X, "/")
    
    if (hyperparameter=="fast"){ 
      DR.sigma <- "median"
    } else {
      DR.sigma <- "auto"
    }
    
    set.seed(seed) 
    DR.Y0A1.Y0A0 <- KLIEP_bw(as.matrix(cbind(Y0s, Xs)[idx1, ]),
                             as.matrix(cbind(Y0s, Xs)[idx0, ]),
                             sigma=DR.sigma,
                             n_y_cols = 1)
    set.seed(seed)
    DR.Yn1A1.Yn1A0 <- KLIEP_bw(as.matrix(cbind(Yn1s, Xs)[idx1, ]),
                               as.matrix(cbind(Yn1s, Xs)[idx0, ]),
                               sigma=DR.sigma,
                               n_y_cols = 1)
    
    ## ----- Construct density ratio -----
    r0 <- matrix(compute_density_ratio_Kernel(as.matrix(cbind(eval_y_vec_s, eval_x_mat_s)),
                                          DR.Y0A1.Y0A0), N)
    rn1 <- matrix(compute_density_ratio_Kernel(as.matrix(cbind(eval_y_vec_s, eval_x_mat_s)),
                                          DR.Yn1A1.Yn1A0), N)
    
    ## ----- Restrict to grid points where density is sufficiently large -----
    
    ## OR_Y0(y,x) = f(y|A=1,X)*f(y_ref|A=0,X) / (f(y_ref|A=1,X)*f(y|A=0,X))
    OR.Y0.alpha0 <- (r0[,1] / r0[,3])
    OR.Yn1.alpha0 <- (r0[,2] / r0[,3])
    
    OR.Y0.alphan1 <- (rn1[,1] / rn1[,3])
    OR.Yn1.alphan1 <- (rn1[,2] / rn1[,3])
    
    range_alpha <- c( OR.Y0.alpha0[idx1] / OR.Y0.alphan1[idx1],
                      OR.Yn1.alpha0[idx1] / OR.Yn1.alphan1[idx1] )
    
    return( stats::quantile(range_alpha, quantile.range) )
    
  } else if (type == "binary"){
    
    ## ----- Binary: use MySL propensity scores (unchanged) -----
    Data.Reset    <- data.frame(cbind(round(A), Yn1, Y0, X))
    colnames(Data.Reset) <- c("A", "Yn1", "Y0", xnames)
    
    set.seed(seed)
    SL.ProbA.givenYn1 <- MySL(Data.Reset, locY = 1, locX = c(2, 3 + 1:d),
                              Ydist = stats::binomial(),
                              SL.list = SL.hpara$SLL, MTRY = SL.hpara$MTRY,
                              MLPL = SL.hpara$MLPL, NMN = SL.hpara$NMN,
                              MLPdecay = SL.hpara$MLPdecay)
    
    set.seed(seed + 1)
    SL.ProbA.givenY0 <- MySL(Data.Reset, locY = 1, locX = c(3, 3 + 1:d),
                             Ydist = stats::binomial(),
                             SL.list = SL.hpara$SLL, MTRY = SL.hpara$MTRY,
                             MLPL = SL.hpara$MLPL, NMN = SL.hpara$NMN,
                             MLPdecay = SL.hpara$MLPdecay)
    
    y_ref     <- 0
    y_grid    <- c(0,1)
    ng <- length(y_grid)
    
    X_rep_g    <- X[rep(seq_len(N), times = ng), , drop = FALSE]
    y_grid_rep <- rep(y_grid, each = N)
    
    Data.Reset.Eval.Y0 <- Data.Reset.Eval.Yn1 <- cbind(y_grid_rep, X_rep_g)
    colnames(Data.Reset.Eval.Y0) <- c("Y0", xnames)
    colnames(Data.Reset.Eval.Yn1) <- c("Yn1", xnames)
    
    Data.Reset.Eval.Y0_ref <- Data.Reset.Eval.Yn1_ref <- cbind(rep(y_ref, N * ng), X_rep_g)
    colnames(Data.Reset.Eval.Y0_ref) <- c("Y0", xnames)
    colnames(Data.Reset.Eval.Yn1_ref) <- c("Yn1", xnames)
    
    r0        <- stats::predict(SL.ProbA.givenY0, newdata = Data.Reset.Eval.Y0, onlySL = TRUE)$pred
    rn1       <- stats::predict(SL.ProbA.givenYn1, newdata = Data.Reset.Eval.Yn1, onlySL = TRUE)$pred
    r0_ref    <- stats::predict(SL.ProbA.givenY0, newdata = Data.Reset.Eval.Y0_ref, onlySL = TRUE)$pred
    rn1_ref   <- stats::predict(SL.ProbA.givenYn1, newdata = Data.Reset.Eval.Yn1_ref, onlySL = TRUE)$pred
    
    ## alpha_0 (Gamma-free)
    alpha0_grid  <- matrix(safe_ratio(r0, r0_ref),N,ng)
    alphan1_grid <- matrix(safe_ratio(rn1, rn1_ref),N,ng)
    range_alpha <- as.numeric(alpha0_grid/alphan1_grid)
    
    return( stats::quantile(range_alpha,quantile.range) )
    
    
  } else {
    stop("type must be 'continuous' or 'binary'")
  }
  
}
















## ============================================================
## PARAMETRIC — NO COVARIATES — SENSITIVITY PARAMETER
## ============================================================

#' @keywords internal
UDID_Parametric_NoX_SensPara <- function(Yn1,
                                         Y0,
                                         A,
                                         type = "continuous",
                                         quantile.range = seq(0,1,by=0.005)) {
  
  N      <- length(Y0)
  idx0   <- which(A == 0)
  idx1   <- which(A == 1)
  pr_A   <- mean(A)
  
  ## Density ratio classifiers (Y only, Gamma-free)
  dr0  <- fit_density_ratio_NoX(Y0[idx1],  Y0[idx0])
  drn1 <- fit_density_ratio_NoX(Yn1[idx1], Yn1[idx0])
  
  if (type == "continuous") {
    
    y_ref  <- stats::median(Y0[idx0])
    
    ## ----- Evaluate precisely at observed Y0 and Yn1 -----
    eval_y_vec <- c(Y0, Yn1, rep(y_ref, N))
    
    r0  <- matrix(predict_density_ratio_NoX(dr0,  eval_y_vec), N)
    rn1 <- matrix(predict_density_ratio_NoX(drn1, eval_y_vec), N)
    
    ## OR_Y0(y) = f(y|A=1)*f(y_ref|A=0) / (f(y_ref|A=1)*f(y|A=0))
    OR.Y0.alpha0  <- r0[,1] / r0[,3]
    OR.Yn1.alpha0 <- r0[,2] / r0[,3]
    
    OR.Y0.alphan1  <- rn1[,1] / rn1[,3]
    OR.Yn1.alphan1 <- rn1[,2] / rn1[,3]
    
    range_alpha <- c( OR.Y0.alpha0[idx1]  / OR.Y0.alphan1[idx1],
                      OR.Yn1.alpha0[idx1] / OR.Yn1.alphan1[idx1] )
    
    return(stats::quantile(range_alpha, quantile.range))
    
  } else if (type == "binary") {
    
    y_ref  <- 0
    y_grid <- c(0, 1)
    ng     <- length(y_grid)
    
    r0      <- predict_density_ratio_NoX(dr0,  y_grid)
    rn1     <- predict_density_ratio_NoX(drn1, y_grid)
    r0_ref  <- predict_density_ratio_NoX(dr0,  rep(y_ref, ng))
    rn1_ref <- predict_density_ratio_NoX(drn1, rep(y_ref, ng))
    
    alpha0_grid  <- matrix(safe_ratio(r0,  r0_ref),  N, ng, byrow = TRUE)
    alphan1_grid <- matrix(safe_ratio(rn1, rn1_ref), N, ng, byrow = TRUE)
    range_alpha <- as.numeric(alpha0_grid/alphan1_grid)
    
    return(stats::quantile(range_alpha, quantile.range))
    
  } else if (type == "poisson") {
    
    y_ref  <- stats::median(Y0[idx0])
    y_grid <- c(0, max(c(Y0[idx1], Yn1[idx1])))
    ng     <- length(y_grid)
    
    r0      <- predict_density_ratio_NoX(dr0,  y_grid)
    rn1     <- predict_density_ratio_NoX(drn1, y_grid)
    r0_ref  <- predict_density_ratio_NoX(dr0,  rep(y_ref, ng))
    rn1_ref <- predict_density_ratio_NoX(drn1, rep(y_ref, ng))
    
    alpha0_grid  <- matrix(safe_ratio(r0,  r0_ref),  N, ng, byrow = TRUE)
    alphan1_grid <- matrix(safe_ratio(rn1, rn1_ref), N, ng, byrow = TRUE)
    range_alpha <- as.numeric(alpha0_grid/alphan1_grid)
    
    return(stats::quantile(range_alpha, quantile.range))
    
  } else {
    stop("type must be 'continuous', 'binary', or 'poisson'")
  }
}


## ============================================================
## NONPARAMETRIC — NO COVARIATES — SENSITIVITY PARAMETER
## ============================================================

UDID_Nonparametric_NoX_SensPara <- function(Yn1,
                                            Y0,
                                            A,
                                            type     = "continuous",
                                            quantile.range = seq(0,1,by=0.005),
                                            seed     = 42,
                                            hyperparameter = "fast",
                                            SL.hpara = list(SLL      = c(1),
                                                            MLPL     = c(2, 4),
                                                            MTRY     = c(1, 2, 3),
                                                            NMN      = 25,
                                                            MLPdecay = c(0.1, 0.001))) {
  
  N      <- length(Y0)
  idx0   <- which(A == 0)
  idx1   <- which(A == 1)
  
  if (type == "continuous") {
    
    ## ----- Standardise Y; use Y0 mean/sd -----
    mu_Y <- mean(Y0);  sd_Y <- stats::sd(Y0)
    Y0s  <- (Y0  - mu_Y) / sd_Y
    Yn1s <- (Yn1 - mu_Y) / sd_Y
    
    y_ref <- stats::median(Y0[idx0])
    
    ## ----- Grid setup -----
    Num.Y.Grid.Basis <- 501
    Y.Grid.Basis <- seq(
      min(c(Y0, Yn1)) - max(stats::sd(Y0), stats::sd(Yn1)),
      max(c(Y0, Yn1) + max(stats::sd(Y0), stats::sd(Yn1))),
      length = Num.Y.Grid.Basis
    )
    dY        <- diff(Y.Grid.Basis)[1]
    grid_cols <- 4:(3 + Num.Y.Grid.Basis)
    
    eval_y_vec <- c(Y0, Yn1,
                    rep(y_ref, N),
                    rep(Y.Grid.Basis, each = N)) 
    
    ## ----- Density ratios via KLIEP_bw -----
    eval_y_vec_s <- (eval_y_vec - mu_Y) / sd_Y
    
    if (hyperparameter == "fast") {
      DR.sigma <- "median"
    } else {
      DR.sigma <- "auto"
    }
    
    set.seed(seed)
    DR.Y0A1.Y0A0 <- KLIEP_bw(matrix(Y0s[idx1]), matrix(Y0s[idx0]),
                             sigma = DR.sigma, n_y_cols = 1)
    
    set.seed(seed)
    DR.Yn1A1.Yn1A0 <- KLIEP_bw(matrix(Yn1s[idx1]), matrix(Yn1s[idx0]),
                               sigma = DR.sigma, n_y_cols = 1)
    
    ## ----- Construct density ratio -----
    r0 <- matrix(compute_density_ratio_Kernel(matrix(eval_y_vec_s),
                                              DR.Y0A1.Y0A0), N)
    rn1 <- matrix(compute_density_ratio_Kernel(matrix(eval_y_vec_s),
                                               DR.Yn1A1.Yn1A0), N)
    
    ## ----- Restrict to grid points where density is sufficiently large -----
    
    ## OR_Y0(y,x) = f(y|A=1,X)*f(y_ref|A=0,X) / (f(y_ref|A=1,X)*f(y|A=0,X))
    OR.Y0.alpha0 <- (r0[,1] / r0[,3])
    OR.Yn1.alpha0 <- (r0[,2] / r0[,3])
    
    OR.Y0.alphan1 <- (rn1[,1] / rn1[,3])
    OR.Yn1.alphan1 <- (rn1[,2] / rn1[,3])
    
    range_alpha <- c( OR.Y0.alpha0[idx1] / OR.Y0.alphan1[idx1],
                      OR.Yn1.alpha0[idx1] / OR.Yn1.alphan1[idx1] )
    
    return( stats::quantile(range_alpha, quantile.range) )
    
  } else if (type == "binary") {
    
    ## Build data with Y only (no X columns)
    Data.Reset <- data.frame(A = round(A), Yn1 = Yn1, Y0 = Y0)
    
    set.seed(seed)
    SL.ProbA.givenYn1 <- MySL(Data.Reset, locY = 1, locX = 2,
                              Ydist = stats::binomial(),
                              SL.list = SL.hpara$SLL, MTRY = SL.hpara$MTRY,
                              MLPL = SL.hpara$MLPL, NMN = SL.hpara$NMN,
                              MLPdecay = SL.hpara$MLPdecay)
    
    set.seed(seed + 1)
    SL.ProbA.givenY0 <- MySL(Data.Reset, locY = 1, locX = 3,
                             Ydist = stats::binomial(),
                             SL.list = SL.hpara$SLL, MTRY = SL.hpara$MTRY,
                             MLPL = SL.hpara$MLPL, NMN = SL.hpara$NMN,
                             MLPdecay = SL.hpara$MLPdecay)
    
    y_ref  <- 0
    y_grid <- c(0, 1)
    ng     <- length(y_grid)
    
    y_grid_rep <- rep(y_grid, each = N)
    
    Data.Reset.Eval.Y0      <- data.frame(Y0  = y_grid_rep)
    Data.Reset.Eval.Yn1     <- data.frame(Yn1 = y_grid_rep)
    Data.Reset.Eval.Y0_ref  <- data.frame(Y0  = rep(y_ref, N * ng))
    Data.Reset.Eval.Yn1_ref <- data.frame(Yn1 = rep(y_ref, N * ng))
    
    r0      <- stats::predict(SL.ProbA.givenY0,  newdata = Data.Reset.Eval.Y0,      onlySL = TRUE)$pred
    rn1     <- stats::predict(SL.ProbA.givenYn1, newdata = Data.Reset.Eval.Yn1,     onlySL = TRUE)$pred
    r0_ref  <- stats::predict(SL.ProbA.givenY0,  newdata = Data.Reset.Eval.Y0_ref,  onlySL = TRUE)$pred
    rn1_ref <- stats::predict(SL.ProbA.givenYn1, newdata = Data.Reset.Eval.Yn1_ref, onlySL = TRUE)$pred
    
    alpha0_grid  <- matrix(safe_ratio(r0,  r0_ref),  N, ng)
    alphan1_grid <- matrix(safe_ratio(rn1, rn1_ref), N, ng)
    range_alpha  <- as.vector(alpha0_grid / alphan1_grid)
    
    return(stats::quantile(range_alpha, quantile.range))
    
  } else {
    stop("type must be 'continuous' or 'binary'")
  }
}