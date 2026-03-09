#' Sensitivity Parameter from Parametric UDID Model (Experimental/Beta Version)
#'
#' Estimates the range of the sensitivity parameter using observed
#' data prior to treatment adoption, based on parametric methods
#' under the UDID framework (Park and Tchetgen Tchetgen, 2026+).
#'
#' @param Yn1 Numeric vector of negative-one-period outcomes (\eqn{Y_{t=-1}}).
#' @param Y0 Numeric vector of pre-treatment outcomes (\eqn{Y_{t=0}}).
#' @param A Binary treatment indicator (0/1).
#' @param X Optional covariate matrix. If \code{NULL}, a no-covariate version is used.
#' @param type Outcome type: \code{"continuous"}, \code{"binary"}, or \code{"poisson"}.
#' @param quantile.range Numeric vector of probabilities in (0,1) specifying
#'   which quantiles to report of the ratio
#'   \eqn{\alpha_{t=0}(y,x) / \alpha_{t=-1}(y,x)}, evaluated at
#'   \eqn{y \in \{Y_{t=0}, Y_{t=-1}\}} among treated units (\eqn{A=1})
#'   and observed \eqn{x}. Defaults to \code{c(0.025, 0.975)},
#'   corresponding to the 2.5th and 97.5th percentiles, giving a 95\%
#'   range that summarizes the plausible spread of the sensitivity
#'   parameter ratio across outcome values and covariate values.
#'
#' @details
#' \code{UDID_Parametric_SensPara} uses an additional pre-pre-treatment period
#' (\eqn{Y_{t=-1}}) to empirically calibrate the sensitivity parameter
#' \eqn{\Gamma}. Under the odds ratio equi-confounding (OREC) assumption,
#' \eqn{\alpha_0(y,x) = \alpha_{-1}(y,x)} for all \eqn{y} and \eqn{x},
#' implying
#' \deqn{
#'   \frac{\alpha_0(y,x)}{\alpha_{-1}(y,x)} = 1
#' }
#' for all \eqn{y,x}. In practice, the OREC assumption may only hold
#' approximately; the distribution of this ratio across outcome and
#' covariate values indicates the plausible range of \eqn{\Gamma}.
#'
#' The density ratios \eqn{f(Y \mid A=1, X)/f(Y \mid A=0, X)} at periods
#' \eqn{t = 0} and \eqn{t = -1} are estimated using logistic regression
#' density ratio classifiers, and the odds ratio functions
#' \eqn{\alpha_0(y,x)} and \eqn{\alpha_{-1}(y,x)} are constructed from
#' these ratios. The function returns selected quantiles of the ratio
#' \eqn{\alpha_0 / \alpha_{-1}} evaluated among treated units.
#'
#' @return A named numeric vector of quantiles of the sensitivity parameter ratio.
#'
#' @references
#' \itemize{
#'   \item Park, C., & Tchetgen Tchetgen, E. (2026+).
#'     A Universal Nonparametric Framework for Difference-in-Differences Analyses.
#'     \url{https://arxiv.org/abs/2212.13641}.
#' }
#'
#' @seealso \code{\link{UDID_Nonparametric_SensPara}} for the nonparametric
#'   version, \code{\link{UDID_Parametric}} for the main estimator.
#'
#' @export
UDID_Parametric_SensPara <- function(Yn1,
                                     Y0,
                                     A,
                                     X    = NULL,
                                     type = "continuous",
                                     quantile.range = c(0.025,0.975)) {
  
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

#' Sensitivity Parameter from Nonparametric UDID Model
#'
#' Estimates the range of the sensitivity parameter using observed
#' data prior to treatment adoption, based on nonparametric methods
#' under the UDID framework (Park and Tchetgen Tchetgen, 2026+).
#'
#' @param Yn1 Numeric vector of negative-one-period outcomes (\eqn{Y_{t=-1}}).
#' @param Y0 Numeric vector of pre-treatment outcomes (\eqn{Y_{t=0}}).
#' @param A Binary treatment indicator (0/1).
#' @param X Optional covariate matrix. If \code{NULL}, a no-covariate version is used.
#' @param type Outcome type: \code{"continuous"} or \code{"binary"}.
#' @param quantile.range Numeric vector of probabilities in (0,1) specifying
#'   which quantiles to report of the ratio
#'   \eqn{\alpha_{t=0}(y,x) / \alpha_{t=-1}(y,x)}, evaluated at
#'   \eqn{y \in \{Y_{t=0}, Y_{t=-1}\}} among treated units (\eqn{A=1})
#'   and observed \eqn{x}. Defaults to \code{c(0.025, 0.975)},
#'   corresponding to the 2.5th and 97.5th percentiles, giving a 95\%
#'   range that summarizes the plausible spread of the sensitivity
#'   parameter ratio across outcome values and covariate values.
#' @param seed Random seed for reproducibility.
#' @param hyperparameter Hyperparameter tuning method: \code{"fast"} or \code{"slow"}.
#'   If \code{"fast"}, the bandwidth of the Gaussian kernel in the density ratio
#'   estimator is determined by the median heuristic (Garreau et al., 2017).
#'   If \code{"slow"}, the bandwidth is selected via cross-validation.
#' @param SL.list Integer vector selecting which learner groups to include (1--9)
#'   in the Super Learner ensemble passed to \code{MySL}.
#'   Available groups: 1 = GLM, 2 = lasso/ridge, 3 = earth,
#'   4 = GAM, 5 = xgboost, 6 = polynomial spline, 7 = random forest,
#'   8 = gbm, 9 = 1-layer MLP.
#'
#' @details
#' \code{UDID_Nonparametric_SensPara} uses an additional pre-pre-treatment
#' period (\eqn{Y_{t=-1}}) to empirically calibrate the sensitivity parameter
#' \eqn{\Gamma}. Under the odds ratio equi-confounding (OREC) assumption,
#' \eqn{\alpha_0(y,x) = \alpha_{-1}(y,x)} for all \eqn{y} and \eqn{x},
#' implying
#' \deqn{
#'   \frac{\alpha_0(y,x)}{\alpha_{-1}(y,x)} = 1
#' }
#' for all \eqn{y,x}. In practice, the OREC assumption may only hold
#' approximately; the distribution of this ratio across outcome and
#' covariate values indicates the plausible range of \eqn{\Gamma}.
#'
#' When the outcome is \strong{continuous}, the density ratios
#' \eqn{f(Y \mid A=1, X)/f(Y \mid A=0, X)} at periods \eqn{t = 0} and
#' \eqn{t = -1} are estimated via the Kullback-Leibler Importance Estimation
#' Procedure (KLIEP; Sugiyama et al., 2007; Nguyen et al., 2007), implemented
#' via the \code{densratio} R package (Makiyama, 2019).
#'
#' When the outcome is \strong{binary}, the density ratios are estimated via
#' an ensemble of machine learning approaches using the Super Learner
#' algorithm (van der Laan et al., 2007), implemented via the
#' \code{SuperLearner} R package (Polley et al., 2025).
#'
#' The odds ratio functions \eqn{\alpha_0(y,x)} and \eqn{\alpha_{-1}(y,x)}
#' are constructed from these ratios, and the function returns selected
#' quantiles of \eqn{\alpha_0 / \alpha_{-1}} evaluated among treated units.
#'
#' @return A named numeric vector of quantiles of the sensitivity parameter ratio.
#'
#' @references
#' \itemize{
#'   \item Park, C., & Tchetgen Tchetgen, E. (2026+).
#'     A Universal Nonparametric Framework for Difference-in-Differences Analyses.
#'     \url{https://arxiv.org/abs/2212.13641}.
#'   \item Garreau, D., Jitkrittum, W., & Kanagawa, M. (2017).
#'     Large sample analysis of the median heuristic.
#'     \url{https://arxiv.org/abs/1707.07269}.
#'   \item Makiyama, K. (2019).
#'     densratio: Density Ratio Estimation. R package version 0.2.1.
#'   \item Nguyen, X., Wainwright, M. J., & Jordan, M. (2007).
#'     Estimating divergence functionals and the likelihood ratio by penalized
#'     convex risk minimization.
#'     \emph{Advances in Neural Information Processing Systems}, 20.
#'   \item Polley, E., LeDell, E., Kennedy, C., Lendle, S., & van der Laan, M. (2025).
#'     SuperLearner: Super Learner Prediction. R package version 2.0-40.
#'   \item Sugiyama, M., Nakajima, S., Kashima, H., Buenau, P., & Kawanabe, M. (2007).
#'     Direct importance estimation with model selection and its application to
#'     covariate shift adaptation.
#'     \emph{Advances in Neural Information Processing Systems}, 20.
#'   \item van der Laan, M. J., Polley, E. C., & Hubbard, A. E. (2007).
#'     Super learner.
#'     \emph{Statistical Applications in Genetics and Molecular Biology}, 6(1).
#' }
#'
#' @seealso \code{\link{UDID_Parametric_SensPara}} for the parametric version,
#'   \code{\link{UDID_Nonparametric}} for the main estimator.
#'
#' @export
UDID_Nonparametric_SensPara <- function(Yn1,
                                        Y0,
                                        A,
                                        X         = NULL,
                                        type      = "continuous",
                                        quantile.range = c(0.025,0.975),
                                        seed      = 42,
                                        hyperparameter = "fast",
                                        SL.list = c(1)) {
  
  ## Dispatch to no-covariate version if X is NULL
  if (is.null(X)) {
    return(UDID_Nonparametric_NoX_SensPara(Yn1 = Yn1, Y0 = Y0, A = A,
                                           type = type, seed = seed,
                                           hyperparameter = hyperparameter,
                                           quantile.range = quantile.range,
                                           SL.list = SL.list))
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
                              SL.list = SL.list)
    
    set.seed(seed + 1)
    SL.ProbA.givenY0 <- MySL(Data.Reset, locY = 1, locX = c(3, 3 + 1:d),
                             Ydist = stats::binomial(),
                             SL.list = SL.list)
    
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
                                         quantile.range = c(0.025,0.975)) {
  
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
                                            quantile.range = c(0.025,0.975),
                                            seed     = 42,
                                            hyperparameter = "fast",
                                            SL.list = c(1)) {
  
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
    return(UDID_Parametric_NoX_SensPara(Yn1 = Yn1, Y0 = Y0, A = A,
                                        type = type, quantile.range = quantile.range))
  } else {
    stop("type must be 'continuous' or 'binary'")
  }
}