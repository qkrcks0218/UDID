#' UDID Nonparametric Estimator
#'
#' Nonparametric universal difference-in-differences estimator with sensitivity
#' analysis (Park & Tchetgen Tchetgen, 2026+). Nuisance models are fitted via 2-fold cross-fitting.
#' ATT is computed for every Gamma value cheaply.
#'
#' @param Y0 Numeric vector of pre-treatment outcomes.
#' @param Y1 Numeric vector of post-treatment outcomes.
#' @param A Binary treatment indicator (0/1).
#' @param X Optional covariate matrix. If \code{NULL}, a no-covariate version is used.
#' @param type Outcome type: \code{"continuous"} or \code{"binary"}.
#' @param log_Gamma_seq Numeric scalar or vector of log(Gamma) sensitivity values (default 0).
#' @param seed Random seed for reproducibility.
#' @param hyperparameter Hyperparameter tuning method: \code{"fast"} or \code{"slow"}.
#'   If \code{"fast"}, the bandwidth of the kernel density estimator is chosen
#'   using the rule-of-thumb method (Silverman, 1986),
#'   and the bandwidth of the Gaussian kernel in the density ratio estimator is
#'   determined by the median heuristic (Garreau et al., 2017).
#'   If \code{"slow"}, these bandwidth parameters are selected via cross-validation.
#' @param SL.list Integer vector selecting which learner groups to include (1--9)
#'   in the Super Learner ensemble passed to \code{MySL}.
#'   Available groups: 1 = GLM, 2 = lasso/ridge, 3 = earth,
#'   4 = GAM, 5 = xgboost, 6 = polynomial spline, 7 = random forest,
#'   8 = gbm, 9 = 1-layer MLP.
#' @param density.report Logical; if \code{TRUE}, return estimated conditional densities.
#' @param hyperparameter.report Logical; if \code{TRUE}, return \code{SL_Library}
#'   (Super Learner ensemble weights) and, for continuous outcomes,
#'   \code{selected_bws} (selected bandwidth parameters). Default is \code{FALSE}.
#' @param verbose Logical; print progress messages.
#'
#' @details
#' \code{UDID_Nonparametric} implements the nonparametric universal
#' difference-in-differences method with 2-fold cross-fitting.
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
#'   \item \eqn{f(Y_0 \mid A=0, X)} is estimated via the nonparametric kernel density
#'   estimator (Hall et al., 2004), implemented via the \code{np} R package
#'   (Hayfield and Racine, 2008).
#'   \item \eqn{f(Y_0 \mid A=1, X)/f(Y_0 \mid A=0, X)} and
#'   \eqn{f(Y_1 \mid A=0, X)/f(Y_0 \mid A=0, X)} are estimated via the
#'   Kullback-Leibler Importance Estimation Procedure (KLIEP;
#'   Sugiyama et al., 2007; Nguyen et al., 2007), implemented via the
#'   \code{densratio} R package (Makiyama, 2019).
#'   \item \eqn{f(A \mid X)} is estimated via an ensemble of machine learning
#'   approaches using the Super Learner algorithm (van der Laan et al., 2007),
#'   implemented via the \code{SuperLearner} R package (Polley et al., 2025).
#' }
#'
#' When the outcome is \strong{binary}, all conditional distributions
#' (\eqn{f(A \mid X)}, \eqn{f(Y_0 \mid A, X)}, and \eqn{f(Y_1 \mid A=0, X)})
#' are estimated via an ensemble of machine learning approaches using the Super
#' Learner algorithm (van der Laan et al., 2007), implemented via the
#' \code{SuperLearner} R package (Polley et al., 2025).
#'
#' Given these nuisance estimates, the average treatment effect on the treated
#' (ATT) is obtained via the efficient influence function.
#'
#' For sensitivity analysis, given the sensitivity parameter \eqn{\Gamma \geq 1},
#' we allow
#' \deqn{
#'   \alpha_1(y,x) \in
#'   \left[\Gamma^{-1} \cdot \alpha_0(y,x),\; \Gamma \cdot \alpha_0(y,x)\right].
#' }
#' When the outcome is \strong{continuous}, we implement the following two-step
#' approach in order to maximize the deviation.
#' \itemize{
#'   \item For each \eqn{X}, obtain
#'   \eqn{y_c(X) := E[Y_1^{(0)} \mid A=1, X]}, which equals
#'   \deqn{
#'     y_c(X) = \frac{E[Y_0\,\alpha_0(Y_0,X) \mid A=1,X]}{E[\alpha_0(Y_0,X) \mid A=1,X]}
#'   }
#'   under the OREC assumption.
#'   \item Then, \eqn{\tilde{\alpha}_1^{UB}(y,X) := \Gamma^{UB}(y)\,\alpha_0(y,X)}
#'   and \eqn{\tilde{\alpha}_1^{LB}(y,X) := \Gamma^{LB}(y)\,\alpha_0(y,X)}, where
#'   \deqn{
#'     \Gamma^{UB}(y) =
#'     \left\{
#'     \begin{array}{ll}
#'       \Gamma^{-1/2} & \text{if } y <    y_c(X) \\
#'       \Gamma^{1/2}  & \text{if } y \geq y_c(X)
#'     \end{array}
#'     \right.
#'     \,, \quad
#'     \Gamma^{LB}(y) =
#'     \left\{
#'     \begin{array}{ll}
#'       \Gamma^{1/2}  & \text{if } y <    y_c(X) \\
#'       \Gamma^{-1/2} & \text{if } y \geq y_c(X)
#'     \end{array}
#'     \right.
#'   }
#'   To satisfy the boundary condition of the odds ratio, i.e.,
#'   \eqn{\alpha_t(y_R, X) = 1}, \eqn{\tilde{\alpha}_1^{UB}} and
#'   \eqn{\tilde{\alpha}_1^{LB}} are normalized as
#'   \deqn{
#'     \alpha_1^{UB}(y,X) =
#'       \frac{\tilde{\alpha}_1^{UB}(y,X)}{\tilde{\alpha}_1^{UB}(y_R,X)},
#'     \qquad
#'     \alpha_1^{LB}(y,X) =
#'       \frac{\tilde{\alpha}_1^{LB}(y,X)}{\tilde{\alpha}_1^{LB}(y_R,X)}.
#'   }
#' }
#'
#' When the outcome is \strong{binary}, the reference value is fixed to
#' \eqn{y_R = 0}. Therefore,
#' \deqn{
#'   \alpha_1^{UB}(y,X) =
#'   \left\{
#'   \begin{array}{ll}
#'     1                           & \text{if } y = 0 \\
#'     \Gamma \cdot \alpha_0(1,X) & \text{if } y = 1
#'   \end{array}
#'   \right.
#'   \,, \quad
#'   \alpha_1^{LB}(y,X) =
#'   \left\{
#'   \begin{array}{ll}
#'     1                               & \text{if } y = 0 \\
#'     \Gamma^{-1} \cdot \alpha_0(1,X) & \text{if } y = 1
#'   \end{array}
#'   \right.
#' }
#'
#' Based on \eqn{\alpha_1^{UB}(y,X)} and \eqn{\alpha_1^{LB}(y,X)}, we obtain
#' the upper and lower bounds on the ATT and their standard errors.
#'
#' @return A named list with the following components:
#'   \describe{
#'     \item{\code{Effect}}{A single-row data frame reporting the estimated
#'       average treatment effect on the treated (ATT) under the UDID method
#'       at \code{log_Gamma = 0}, its asymptotic standard error (\code{SE}),
#'       and its multiplier bootstrap standard error (\code{Boot.SE}).}
#'     \item{\code{EIF}}{A data frame with one row per observation containing
#'       the uncentered efficient influence function (\code{Uncentered EIF})
#'       and the centered efficient influence function (\code{Centered EIF}),
#'       both evaluated at \code{Gamma = 1} (i.e., \code{log_Gamma = 0}).}
#'     \item{\code{Sensitivity}}{A data frame with one row per entry of
#'       \code{log_Gamma_seq}, reporting the lower and upper sensitivity bounds
#'       on the ATT (\code{ATT_LB}, \code{ATT_UB}) along with their asymptotic
#'       and multiplier bootstrap standard errors (\code{SE_LB}, \code{SE_UB},
#'       \code{SE.Boot_LB}, \code{SE.Boot_UB}).}
#'     \item{\code{SL_Library}}{(only when \code{hyperparameter.report = TRUE})
#'       A data frame reporting the Super Learner ensemble weights assigned to
#'       each candidate algorithm in each cross-fitting fold. For continuous
#'       outcomes the columns cover \code{Pr(A|X)}; for binary outcomes they
#'       additionally cover \code{Pr(Y0|A=0,X)}, \code{Pr(Y0|A=1,X)}, and
#'       \code{Pr(Y1|A=0,X)}.}
#'     \item{\code{selected_bws}}{(only when \code{hyperparameter.report = TRUE}
#'       and \code{type = "continuous"}) A data frame reporting the selected
#'       bandwidth parameters, with one row for the outcome variable Y and one
#'       row per covariate X1, ..., Xd. Columns cover the kernel density
#'       estimate \eqn{f(Y_0|A=0,X)} (separate Y and X bandwidths from
#'       \code{npcdens}) and the two KLIEP density ratio estimates
#'       \eqn{f(Y_0|A=1,X)/f(Y_0|A=0,X)} and
#'       \eqn{f(Y_1|A=0,X)/f(Y_0|A=0,X)} (a single Gaussian kernel bandwidth
#'       \eqn{\sigma} shared across all dimensions), each reported for the two
#'       cross-fitting folds.}
#'     \item{\code{Density}}{(only when \code{density.report = TRUE}) A list of
#'       matrices \code{fY0A0}, \code{fY0A1}, \code{fY1A0}, \code{fY1A1} and
#'       a vector \code{Y.Grid}. Each matrix has one row per observation and
#'       one column per grid point of Y (501 equally spaced points for
#'       continuous outcomes; two points \{0, 1\} for binary outcomes).
#'       Matrix entries give the estimated conditional density
#'       \eqn{f(y \mid A = a, X_i)} at each grid point. For binary outcomes,
#'       \code{fY1A1} is derived from the UDID parallel-trends odds-ratio
#'       assumption at \eqn{\Gamma = 1}.}
#'   }
#'
#' @references
#' \itemize{
#'   \item Park, C., & Tchetgen Tchetgen, E. (2026+).
#'     A Universal Nonparametric Framework for Difference-in-Differences Analyses.
#'     \url{https://arxiv.org/abs/2212.13641}.
#'   \item Silverman, B. W. (1986).
#'     \emph{Density Estimation for Statistics and Data Analysis}.
#'     London: Chapman and Hall.
#'   \item Garreau, D., Jitkrittum, W., & Kanagawa, M. (2017).
#'     Large sample analysis of the median heuristic.
#'     \url{https://arxiv.org/abs/1707.07269}.
#'   \item Hall, P., Racine, J. S., & Li, Q. (2004).
#'     Cross-validation and the estimation of conditional probability densities.
#'     \emph{Journal of the American Statistical Association}, 99, 1015--1026.
#'   \item Sugiyama, M., Nakajima, S., Kashima, H., Buenau, P., & Kawanabe, M. (2007).
#'     Direct importance estimation with model selection and its application to
#'     covariate shift adaptation.
#'     \emph{Advances in Neural Information Processing Systems}, 20.
#'   \item Nguyen, X., Wainwright, M. J., & Jordan, M. (2007).
#'     Estimating divergence functionals and the likelihood ratio by penalized
#'     convex risk minimization.
#'     \emph{Advances in Neural Information Processing Systems}, 20.
#'   \item van der Laan, M. J., Polley, E. C., & Hubbard, A. E. (2007).
#'     Super learner.
#'     \emph{Statistical Applications in Genetics and Molecular Biology}, 6(1).
#'   \item Hayfield, T., & Racine, J. S. (2008).
#'     Nonparametric econometrics: The np package.
#'     \emph{Journal of Statistical Software}, 27, 1--32.
#'   \item Makiyama, K. (2019).
#'     densratio: Density Ratio Estimation. R package version 0.2.1.
#'   \item Polley, E., LeDell, E., Kennedy, C., Lendle, S., & van der Laan, M. (2025).
#'     SuperLearner: Super Learner Prediction. R package version 2.0-40.
#' }
#'
#' @export

UDID_Nonparametric <- function(Y0,
                               Y1,
                               A,
                               X                     = NULL,
                               type                  = "continuous",
                               log_Gamma_seq         = 0,
                               seed                  = 42,
                               hyperparameter        = "fast",
                               SL.list               = c(1),
                               density.report        = FALSE,
                               hyperparameter.report = FALSE,
                               verbose               = FALSE) {

  if (is.null(X)) {
    return(UDID_Nonparametric_NoX(Y0 = Y0, Y1 = Y1, A = A,
                                  type = type, log_Gamma_seq = log_Gamma_seq,
                                  seed = seed, hyperparameter = hyperparameter,
                                  SL.list = SL.list, density.report = density.report,
                                  verbose = verbose))
  }

  Gamma_seq <- exp(log_Gamma_seq)
  N  <- length(Y0)
  X  <- as.matrix(X)
  nG <- length(Gamma_seq)

  set.seed(seed)
  idx1     <- sample(which(A == 1), floor(sum(A == 1) / 2))
  idx0     <- sample(which(A == 0), ceiling(sum(A == 0) / 2))
  SS.Index <- list(c(idx1, idx0), setdiff(seq_len(N), c(idx1, idx0)))

  uncEIF_UB_full   <- replicate(nG, numeric(N), simplify = FALSE)
  uncEIF_LB_full   <- replicate(nG, numeric(N), simplify = FALSE)
  uncEIF_base_full <- numeric(N)

  if (density.report) {
    Num.Grid <- if (type == "continuous") 501L else 2L
    density_full <- list(fY0A0  = matrix(NA_real_, N, Num.Grid),
                         fY0A1  = matrix(NA_real_, N, Num.Grid),
                         fY1A0  = matrix(NA_real_, N, Num.Grid),
                         fY1A1  = matrix(NA_real_, N, Num.Grid),
                         Y.Grid = NULL)
  }

  if (hyperparameter.report) {
    fold_sl_info <- vector("list", 2)
    if (type == "continuous") fold_bw_info <- vector("list", 2)
  }

  for (ss in 1:2) {
    train <- SS.Index[[ss]]; eval <- SS.Index[[3 - ss]]
    fold_results <- .UDID_fold(
      Y0[train], Y1[train], A[train], X[train, , drop = FALSE],
      Y0[eval],  Y1[eval],  A[eval],  X[eval,  , drop = FALSE],
      type, SL.list, Gamma_seq, seed = seed * 10L + ss, fold_id = ss,
      hyperparameter = hyperparameter, density.report = density.report,
      hyperparameter.report = hyperparameter.report, verbose = verbose)
    for (gi in seq_len(nG)) {
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
    if (hyperparameter.report) {
      fold_sl_info[[ss]] <- fold_results$SL_info
      if (type == "continuous") fold_bw_info[[ss]] <- fold_results$bw_info
    }
  }

  pr_A    <- mean(A)
  results <- vector("list", nG)
  base_gi <- which(Gamma_seq == 1)

  for (gi in seq_len(nG)) {
    uncEIF_UB  <- uncEIF_UB_full[[gi]] / pr_A
    tau_UB     <- mean(uncEIF_UB)
    cEIF_UB    <- uncEIF_UB - A * tau_UB / pr_A
    SE_UB      <- stats::sd(cEIF_UB) / sqrt(N)
    set.seed(seed + gi)
    SE.Boot_UB <- mboot_sd_cpp(uncEIF_UB, 1000L)   # BLAS path: avoids 320 MB matrix

    uncEIF_LB  <- uncEIF_LB_full[[gi]] / pr_A
    tau_LB     <- mean(uncEIF_LB)
    cEIF_LB    <- uncEIF_LB - A * tau_LB / pr_A
    SE_LB      <- stats::sd(cEIF_LB) / sqrt(N)
    set.seed(seed + gi + 1000L)
    SE.Boot_LB <- mboot_sd_cpp(uncEIF_LB, 1000L)

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

  eif_unc <- uncEIF_base_full / pr_A
  eif_tau <- mean(eif_unc)
  eif_cen <- eif_unc - A * eif_tau / pr_A
  eif_df  <- data.frame(`Uncentered EIF` = eif_unc, `Centered EIF` = eif_cen,
                        check.names = FALSE)

  RESULT <- list(Effect = effect_row, EIF = eif_df, Sensitivity = sens_df)

  if (hyperparameter.report) {
    sl_i1    <- fold_sl_info[[1]]; sl_i2 <- fold_sl_info[[2]]
    sl_names <- sl_i1$ProbA$names
    sl_table <- data.frame(
      `SuperLearner Library` = sl_names,
      `Fold 1: Pr(A|X)`     = as.numeric(sl_i1$ProbA$coef),
      `Fold 2: Pr(A|X)`     = as.numeric(sl_i2$ProbA$coef),
      check.names = FALSE, stringsAsFactors = FALSE)
    if (type == "binary") {
      sl_table[["Fold 1: Pr(Y0|A=0,X)"]] <- as.numeric(sl_i1$ProbY0A0$coef)
      sl_table[["Fold 2: Pr(Y0|A=0,X)"]] <- as.numeric(sl_i2$ProbY0A0$coef)
      sl_table[["Fold 1: Pr(Y0|A=1,X)"]] <- as.numeric(sl_i1$ProbY0A1$coef)
      sl_table[["Fold 2: Pr(Y0|A=1,X)"]] <- as.numeric(sl_i2$ProbY0A1$coef)
      sl_table[["Fold 1: Pr(Y1|A=0,X)"]] <- as.numeric(sl_i1$ProbY1A0$coef)
      sl_table[["Fold 2: Pr(Y1|A=0,X)"]] <- as.numeric(sl_i2$ProbY1A0$coef)
    }
    RESULT$SL_Library <- sl_table
    if (type == "continuous") {
      bw_i1 <- fold_bw_info[[1]]; bw_i2 <- fold_bw_info[[2]]
      bw_table <- data.frame(
        `Fold 1: f(Y0|A=0,X)`              = c(bw_i1$fY0A0$ybw, as.numeric(bw_i1$fY0A0$xbw)),
        `Fold 2: f(Y0|A=0,X)`              = c(bw_i2$fY0A0$ybw, as.numeric(bw_i2$fY0A0$xbw)),
        `Fold 1: f(Y0|A=1,X)/f(Y0|A=0,X)` = bw_i1$DR_Y0A1_Y0A0$sigma,
        `Fold 2: f(Y0|A=1,X)/f(Y0|A=0,X)` = bw_i2$DR_Y0A1_Y0A0$sigma,
        `Fold 1: f(Y1|A=0,X)/f(Y0|A=0,X)` = bw_i1$DR_Y1A0_Y0A0$sigma,
        `Fold 2: f(Y1|A=0,X)/f(Y0|A=0,X)` = bw_i2$DR_Y1A0_Y0A0$sigma,
        row.names = c("Y", paste0("X", seq_len(ncol(X)))), check.names = FALSE)
      RESULT$selected_bws <- bw_table
    }
  }

  if (density.report) RESULT$Density <- density_full
  return(RESULT)
}


## ---- Fold worker (with covariates) ------------------------------------------

.UDID_fold <- function(Y0, Y1, A, X,
                       Y0.Eval, Y1.Eval, A.Eval, X.Eval,
                       type, SL.list, Gamma_seq, seed, fold_id = NULL,
                       hyperparameter        = "fast",
                       density.report        = FALSE,
                       hyperparameter.report = FALSE,
                       verbose               = FALSE) {

  N.Test           <- length(Y0.Eval)
  Num.Y.Grid.Basis <- 501
  d                <- ncol(X)

  Y.Grid.Basis <- seq(
    min(c(Y0, Y1[A == 0])) - max(stats::sd(Y0), stats::sd(Y1[A == 0])),
    max(c(Y0, Y1[A == 0])) + max(stats::sd(Y0), stats::sd(Y1[A == 0])),
    length = Num.Y.Grid.Basis
  )
  dY        <- diff(Y.Grid.Basis)[1]
  grid_cols <- 4:(3 + Num.Y.Grid.Basis)

  ## ----------------------------------------------------------
  ## Nuisance model estimation (Gamma-free)
  ## ----------------------------------------------------------

  if (type == "continuous") {

    eval_y_vec <- c(Y0.Eval, Y1.Eval,
                    rep(stats::median(Y0[A == 0]), N.Test),
                    rep(Y.Grid.Basis, each = N.Test))
    eval_x_mat <- X.Eval[rep(seq_len(N.Test), times = 3 + Num.Y.Grid.Basis), , drop = FALSE]

    set.seed(seed)
    if (verbose) cat(sprintf("Fold %d: estimating the density f(Y_{t=0}|A=0,X)\n", fold_id))
    if (hyperparameter == "fast") {
      npbw.Y0A0 <- np::npcdensbw(xdat = X[A == 0, ], ydat = Y0[A == 0],
                                 bwmethod = "normal-reference")
    } else {
      npbw.Y0A0 <- np::npcdensbw(xdat = X[A == 0, ], ydat = Y0[A == 0])
    }
    kde_mat      <- matrix(stats::fitted(np::npcdens(npbw.Y0A0,
                                                     exdat = eval_x_mat,
                                                     eydat = eval_y_vec)), N.Test)
    Raw.Y0A0.KDE <- kde_mat / (rowSums(kde_mat[, grid_cols]) * dY)

    DR.sigma <- if (hyperparameter == "fast") "median" else "auto"

    set.seed(seed + 1)
    if (verbose) cat(sprintf("Fold %d: estimating the density ratio f(Y_{t=0}|A=1,X)/f(Y_{t=0}|A=0,X)\n", fold_id))
    mu_Y <- mean(Y0);  sd_Y <- stats::sd(Y0)
    mu_X <- colMeans(X);  sd_X <- apply(X, 2, sd);  sd_X[sd_X == 0] <- 1
    Y0s  <- (Y0 - mu_Y) / sd_Y;  Y1s <- (Y1 - mu_Y) / sd_Y
    Xs   <- sweep(sweep(X, 2, mu_X), 2, sd_X, "/")
    eval_y_vec_s <- (eval_y_vec - mu_Y) / sd_Y
    eval_x_mat_s <- sweep(sweep(eval_x_mat, 2, mu_X), 2, sd_X, "/")

    DR.Y0A1.Y0A0 <- KLIEP_bw(as.matrix(cbind(Y0s, Xs)[A == 1, ]),
                              as.matrix(cbind(Y0s, Xs)[A == 0, ]),
                              n_y_cols = 1, sigma = DR.sigma, verbose = verbose)
    Raw.Y0A1.DR  <- Raw.Y0A0.KDE *
      matrix(compute_density_ratio_Kernel(as.matrix(cbind(eval_y_vec_s, eval_x_mat_s)),
                                          DR.Y0A1.Y0A0), N.Test)
    Raw.Y0A1.DR  <- Raw.Y0A1.DR / (rowSums(Raw.Y0A1.DR[, grid_cols]) * dY)

    set.seed(seed + 2)
    if (verbose) cat(sprintf("Fold %d: estimating the density ratio f(Y_{t=1}|A=0,X)/f(Y_{t=0}|A=0,X)\n", fold_id))
    DR.Y1A0.Y0A0 <- KLIEP_bw(as.matrix(cbind(Y1s, Xs)[A == 0, ]),
                              as.matrix(cbind(Y0s, Xs)[A == 0, ]),
                              n_y_cols = 1, sigma = DR.sigma, verbose = verbose)
    Raw.Y1A0.DR  <- Raw.Y0A0.KDE *
      matrix(compute_density_ratio_Kernel(as.matrix(cbind(eval_y_vec_s, eval_x_mat_s)),
                                          DR.Y1A0.Y0A0), N.Test)
    Raw.Y1A0.DR  <- Raw.Y1A0.DR / (rowSums(Raw.Y1A0.DR[, grid_cols]) * dY)

    col_Y0 <- 1L; col_Y1 <- 2L; col_base <- 3L

    xnms <- sprintf("X%0.5d", 1:d)
    set.seed(seed + 3)
    if (verbose) cat(sprintf("Fold %d: estimating the propensity score Pr(A|X)\n", fold_id))
    ProbA1.Fit <- MySL(stats::setNames(data.frame(cbind(A, X)), c("A", xnms)),
                       locY = 1, locX = 1 + 1:d, Ydist = stats::binomial(),
                       SL.list = SL.list)
    Est.Marginal.ProbA <- stats::predict(ProbA1.Fit,
                                         newdata = stats::setNames(data.frame(X.Eval), xnms),
                                         onlySL  = TRUE)$pred

    if (hyperparameter.report) {
      SL_info <- list(ProbA = list(names = colnames(ProbA1.Fit$library.predict),
                                   coef  = ProbA1.Fit$coef))
      bw_info <- list(fY0A0        = list(ybw = as.numeric(npbw.Y0A0$ybw),
                                          xbw = as.numeric(npbw.Y0A0$xbw)),
                      DR_Y0A1_Y0A0 = list(sigma = DR.Y0A1.Y0A0$kernel_info$sigma),
                      DR_Y1A0_Y0A0 = list(sigma = DR.Y1A0.Y0A0$kernel_info$sigma))
    }

    if (verbose) cat(sprintf("Fold %d: calculating the efficient influence function\n", fold_id))

    ## Gamma-free quantities — vectorised pmin/pmax (no per-element apply/median)
    Joint.Y0A.Y0 <- cbind(
      Raw.Y0A1.DR[, col_Y0]    * Est.Marginal.ProbA,
      Raw.Y0A1.DR[, col_base]  * Est.Marginal.ProbA,
      Raw.Y0A0.KDE[, col_Y0]   * (1 - Est.Marginal.ProbA),
      Raw.Y0A0.KDE[, col_base] * (1 - Est.Marginal.ProbA)
    )
    Joint.Y0A.Y1 <- cbind(
      Raw.Y0A1.DR[, col_Y1]    * Est.Marginal.ProbA,
      Raw.Y0A1.DR[, col_base]  * Est.Marginal.ProbA,
      Raw.Y0A0.KDE[, col_Y1]   * (1 - Est.Marginal.ProbA),
      Raw.Y0A0.KDE[, col_base] * (1 - Est.Marginal.ProbA)
    )
    hat.OR0.alpha0 <- pmin(pmax(Joint.Y0A.Y0[,1] * Joint.Y0A.Y0[,4] /
                                  (Joint.Y0A.Y0[,2] * Joint.Y0A.Y0[,3]), 0), 1e10)
    hat.OR1.alpha0 <- pmin(pmax(Joint.Y0A.Y1[,1] * Joint.Y0A.Y1[,4] /
                                  (Joint.Y0A.Y1[,2] * Joint.Y0A.Y1[,3]), 0), 1e10)
    hat.BetaA0 <- Joint.Y0A.Y0[,2] / Joint.Y0A.Y0[,4]
    DR.base    <- Raw.Y1A0.DR[, col_Y0] / Raw.Y0A0.KDE[, col_Y0]

    OR.Grid.alpha0 <- pmin(pmax(
      Raw.Y0A1.DR[, grid_cols] *
        matrix(Raw.Y0A0.KDE[, col_base], N.Test, Num.Y.Grid.Basis) /
        (matrix(Raw.Y0A1.DR[, col_base], N.Test, Num.Y.Grid.Basis) *
           Raw.Y0A0.KDE[, grid_cols]),
      0), 1e10)
    Cond.Density <- Raw.Y1A0.DR[, grid_cols] * dY

  } else {
    ## ---- Binary outcome branch ----
    xnms       <- sprintf("X%0.5d", 1:d)
    Data.Reset <- data.frame(cbind(round(Y0), round(Y1), round(A), X))
    colnames(Data.Reset) <- c("Y0", "Y1", "A", xnms)
    Data.Reset.A0      <- Data.Reset[Data.Reset$A == 0, ]
    Data.Reset.Eval    <- stats::setNames(data.frame(X.Eval), xnms)
    Data.Reset.A0.Eval <- stats::setNames(data.frame(cbind(0, X.Eval)), c("A", xnms))
    Data.Reset.A1.Eval <- stats::setNames(data.frame(cbind(1, X.Eval)), c("A", xnms))

    set.seed(seed + 4)
    if (verbose) cat(sprintf("Fold %d: estimating the propensity score Pr(A|X)\n", fold_id))
    SL.ProbA <- MySL(Data.Reset, locY = 3, locX = 3 + 1:d,
                     Ydist = stats::binomial(), SL.list = SL.list)
    set.seed(seed + 5)
    if (verbose) cat(sprintf("Fold %d: estimating the outcome regression Pr(Y_{t=0}|A,X)\n", fold_id))
    SL.ProbY0gA <- MySL(Data.Reset, locY = 1, locX = c(3, 3 + 1:d),
                        Ydist = stats::binomial(), SL.list = SL.list)
    set.seed(seed + 7)
    if (verbose) cat(sprintf("Fold %d: estimating the outcome regression Pr(Y_{t=1}|A=0,X)\n", fold_id))
    SL.ProbY1gA0 <- MySL(Data.Reset.A0, locY = 2, locX = 3 + 1:d,
                         Ydist = stats::binomial(), SL.list = SL.list)

    if (hyperparameter.report) {
      SL_info <- list(
        ProbA    = list(names = colnames(SL.ProbA$library.predict),     coef = SL.ProbA$coef),
        ProbY0A0 = list(names = colnames(SL.ProbY0gA$library.predict),  coef = SL.ProbY0gA$coef),
        ProbY0A1 = list(names = colnames(SL.ProbY0gA$library.predict),  coef = SL.ProbY0gA$coef),
        ProbY1A0 = list(names = colnames(SL.ProbY1gA0$library.predict), coef = SL.ProbY1gA0$coef)
      )
    }

    if (verbose) cat(sprintf("Fold %d: calculating the efficient influence function\n", fold_id))
    bound_prob   <- function(x, tol = 1e-4) pmax(pmin(x, 1 - tol), tol)
    Est.PrA1     <- bound_prob(stats::predict(SL.ProbA,     newdata = Data.Reset.Eval,    onlySL = TRUE)$pred)
    Est.PrY01gA0 <- bound_prob(stats::predict(SL.ProbY0gA,  newdata = Data.Reset.A0.Eval, onlySL = TRUE)$pred)
    Est.PrY01gA1 <- bound_prob(stats::predict(SL.ProbY0gA,  newdata = Data.Reset.A1.Eval, onlySL = TRUE)$pred)
    Est.PrY11gA0 <- bound_prob(stats::predict(SL.ProbY1gA0, newdata = Data.Reset.Eval,    onlySL = TRUE)$pred)

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

    if (density.report) {
      bin_fY0A0  <- cbind(1 - Est.PrY01gA0, Est.PrY01gA0)
      bin_fY0A1  <- cbind(1 - Est.PrY01gA1, Est.PrY01gA1)
      bin_fY1A0  <- cbind(1 - Est.PrY11gA0, Est.PrY11gA0)
      denom_Y1A1 <- (1 - Cond.Density) + hat.OR.x.alpha0 * Cond.Density
      bin_fY1A1  <- cbind((1 - Cond.Density) / denom_Y1A1,
                          hat.OR.x.alpha0 * Cond.Density / denom_Y1A1)
    }
  }

  ## ----------------------------------------------------------
  ## Sweep over Gamma_seq — C++-accelerated inner loop
  ## ----------------------------------------------------------

  pr_A.eval <- mean(A.Eval)

  .fold_eif <- function(Gamma, direction) {

    if (type == "continuous") {

      if (Gamma == 1) {
        ## C++: single pass, no temporary N x 501 matrices
        rs             <- rowsums_grid_cpp(OR.Grid.alpha0, Cond.Density, Y.Grid.Basis)
        E_alpha1       <- rs$E_alpha
        E_Y1alpha1     <- rs$E_Yalpha
        hat.OR1        <- hat.OR1.alpha0
        hat.OR0.alpha1 <- hat.OR0.alpha0

      } else {
        ## C++: scale + rowSums fused — saves OR.Grid.alpha1, scale_grid,
        ##      mu_base_grid and two rowSums products (5 temp N x 501 mats)
        rs         <- sens_rowsums_cpp(OR.Grid.alpha0, Cond.Density, Y.Grid.Basis,
                                       mu_base_vec, log(Gamma) / 2, direction == "UB")
        E_alpha1   <- rs$E_alpha
        E_Y1alpha1 <- rs$E_Yalpha
        g_up   <- Gamma ^ 0.5
        g_down <- Gamma ^ (-0.5)
        if (direction == "UB") {
          hat.OR0.alpha1 <- hat.OR0.alpha0 * ifelse(Y0.Eval > mu_base_vec, g_up, g_down)
          hat.OR1        <- hat.OR1.alpha0 * ifelse(Y1.Eval > mu_base_vec, g_up, g_down)
        } else {
          hat.OR0.alpha1 <- hat.OR0.alpha0 * ifelse(Y0.Eval > mu_base_vec, g_down, g_up)
          hat.OR1        <- hat.OR1.alpha0 * ifelse(Y1.Eval > mu_base_vec, g_down, g_up)
        }
      }

      hat.Mu1.plugin    <- E_Y1alpha1 / E_alpha1
      hat.BetaA1.plugin <- Est.Marginal.ProbA / (1 - Est.Marginal.ProbA) / E_alpha1
      hat.DR.plugin     <- (A.Eval / (hat.BetaA0 * hat.OR0.alpha0) + (1 - A.Eval)) *
                           hat.OR0.alpha1 * DR.base * hat.BetaA1.plugin

    } else {
      ## Binary
      if (Gamma == 1) {
        hat.OR.x.alpha1.loc <- hat.OR.x.alpha0; alpha1_at_0 <- 1
      } else {
        if (direction == "UB") {
          s1 <- sens_scale_UB(rep(1, N.Test), mu_base_vec, Gamma)
          s0 <- sens_scale_UB(rep(0, N.Test), mu_base_vec, Gamma)
        } else {
          s1 <- sens_scale_LB(rep(1, N.Test), mu_base_vec, Gamma)
          s0 <- sens_scale_LB(rep(0, N.Test), mu_base_vec, Gamma)
        }
        hat.OR.x.alpha1.loc <- hat.OR.x.alpha0 * s1; alpha1_at_0 <- 1 * s0
      }
      E_alpha1          <- Cond.Density * hat.OR.x.alpha1.loc + (1 - Cond.Density) * alpha1_at_0
      hat.OR1           <- ifelse(Y1.Eval == 1, hat.OR.x.alpha1.loc, alpha1_at_0)
      hat.Mu1.plugin    <- hat.OR.x.alpha1.loc * Cond.Density / E_alpha1
      hat.BetaA1.plugin <- Est.Marginal.ProbA / (1 - Est.Marginal.ProbA) / E_alpha1
      hat.DR.plugin     <- (A.Eval / (beta0.eval * hat.OR.x.alpha0 ^ Y0.Eval) + (1 - A.Eval)) *
                           ifelse(Y0.Eval == 1, hat.OR.x.alpha1.loc, alpha1_at_0) *
                           DR.base * hat.BetaA1.plugin
    }

    uncEIF <- A.Eval * Y1.Eval / pr_A.eval -
      EIF_Continuous(Y1.Eval, Y0.Eval, A.Eval,
                     hat.OR1, hat.BetaA1.plugin, hat.Mu1.plugin,
                     hat.DR.plugin) / pr_A.eval
    list(eif = uncEIF * pr_A.eval, Mu1 = hat.Mu1.plugin)
  }

  ## Baseline at Gamma = 1; mu_base_vec (N-vector) is all sens_rowsums_cpp needs
  base_res    <- .fold_eif(Gamma = 1, direction = "none")
  mu_base_vec <- base_res$Mu1

  results <- vector("list", length(Gamma_seq))
  for (gi in seq_along(Gamma_seq)) {
    Gamma <- Gamma_seq[gi]
    if (Gamma == 1) {
      results[[gi]] <- list(UB = base_res$eif, LB = base_res$eif)
    } else {
      res_UB <- .fold_eif(Gamma, "UB")
      res_LB <- .fold_eif(Gamma, "LB")
      results[[gi]] <- list(UB = res_UB$eif, LB = res_LB$eif)
    }
  }
  results$base_eif <- base_res$eif

  if (density.report && type == "continuous") {
    results$fY0A0  <- Raw.Y0A0.KDE[, grid_cols]
    results$fY0A1  <- Raw.Y0A1.DR[, grid_cols]
    results$fY1A0  <- Raw.Y1A0.DR[, grid_cols]
    fY1A1          <- OR.Grid.alpha0 * Raw.Y1A0.DR[, grid_cols] * dY
    results$fY1A1  <- fY1A1 / matrix(rowSums(fY1A1), N.Test, Num.Y.Grid.Basis)
    results$Y.Grid <- Y.Grid.Basis
  } else if (density.report && type == "binary") {
    results$fY0A0  <- bin_fY0A0; results$fY0A1 <- bin_fY0A1
    results$fY1A0  <- bin_fY1A0; results$fY1A1 <- bin_fY1A1
    results$Y.Grid <- c(0L, 1L)
  }

  if (hyperparameter.report) {
    results$SL_info <- SL_info
    if (type == "continuous") results$bw_info <- bw_info
  }

  results
}
