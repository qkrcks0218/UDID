#' Sensitivity Bounds
#'
#' Compute ATT and confidence interval bounds from a sensitivity analysis
#' output (Park & Tchetgen Tchetgen, 2026+).
#'
#' @param output Output from \code{\link{UDID_Nonparametric}} or
#'   \code{\link{UDID_Parametric}}. Accepts the full estimator output (a list
#'   with a \code{Sensitivity} element), an already-processed bounds object
#'   (with \code{bounds} and \code{crossings}), or a bare \code{Sensitivity}
#'   data frame.
#' @param alpha Confidence level (default 0.95).
#'
#' @details
#' \code{UDID_Sensitivity_Bounds} takes the sensitivity analysis output from
#' \code{\link{UDID_Nonparametric}} or \code{\link{UDID_Parametric}} and
#' computes:
#' \itemize{
#'   \item Lower and upper bounds on the ATT (\code{ATT_LB}, \code{ATT_UB})
#'   as a function of \eqn{\log \Gamma}.
#'   \item Pointwise \eqn{(1-\alpha)} confidence interval bounds
#'   (\code{CI_LB}, \code{CI_UB}), defined as
#'   \deqn{
#'     CI_{LB} = ATT_{LB} - z_{\alpha/2} \cdot SE_{LB}, \qquad
#'     CI_{UB} = ATT_{UB} + z_{\alpha/2} \cdot SE_{UB}.
#'   }
#'   \item Zero-crossing values: the interpolated \eqn{\log \Gamma} at which
#'   each of \code{ATT_LB}, \code{ATT_UB}, \code{CI_LB}, and \code{CI_UB}
#'   crosses zero. These indicate the degree of unmeasured confounding
#'   required to nullify the treatment effect or its statistical significance.
#' }
#'
#' @return A list with two components:
#'   \describe{
#'     \item{\code{bounds}}{A data frame with columns \code{log_Gamma},
#'       \code{ATT_LB}, \code{SE_LB}, \code{ATT_UB}, \code{SE_UB},
#'       \code{CI_LB}, and \code{CI_UB}.}
#'     \item{\code{crossings}}{A named list with the interpolated
#'       \eqn{\log \Gamma} values at which \code{ATT_LB}, \code{ATT_UB},
#'       \code{CI_LB}, and \code{CI_UB} cross zero (\code{NA} if no
#'       crossing is found).}
#'   }
#'
#' @references
#' \itemize{
#'   \item Park, C., & Tchetgen Tchetgen, E. (2026+).
#'     A Universal Nonparametric Framework for Difference-in-Differences Analyses.
#'     \url{https://arxiv.org/abs/2212.13641}.
#' }
#'
#' @seealso \code{\link{UDID_Sensitivity_Plot}} to visualize these bounds,
#'   \code{\link{UDID_Nonparametric}} and \code{\link{UDID_Parametric}} for
#'   the estimators that produce the input.
#'
#' @export
UDID_Sensitivity_Bounds <- function(output, alpha = 0.95) {
  ## Accept: (1) full estimator output with $Sensitivity,
  ##         (2) already-processed bounds object with $bounds,
  ##         (3) a bare Sensitivity data.frame
  if (!is.null(output$bounds) && !is.null(output$crossings)) {
    ## Already processed — return as-is
    return(output)
  }
  if (!is.null(output$Sensitivity)) output <- output$Sensitivity
  
  z  <- stats::qnorm(1 - (1 - alpha) / 2)
  bounds <- data.frame(
    log_Gamma = output$log_Gamma,
    ATT_LB    = output$ATT_LB,
    SE_LB     = output$SE_LB,
    ATT_UB    = output$ATT_UB,
    SE_UB     = output$SE_UB,
    CI_LB     = output$ATT_LB - z * output$SE_LB,
    CI_UB     = output$ATT_UB + z * output$SE_UB
  )
  
  cross <- function(v) {
    x <- bounds$log_Gamma
    if (length(v) < 2) return(NA_real_)
    for (i in seq_len(length(v) - 1)) {
      if (!is.na(v[i]) && !is.na(v[i+1]) &&
          ((v[i] >= 0 && v[i+1] < 0) || (v[i] <= 0 && v[i+1] > 0)))
        return(x[i] + (0 - v[i]) * (x[i+1] - x[i]) / (v[i+1] - v[i]))
    }
    NA_real_
  }
  
  list(
    bounds    = bounds,
    crossings = list(
      ATT_LB = cross(bounds$ATT_LB),
      ATT_UB = cross(bounds$ATT_UB),
      CI_LB  = cross(bounds$CI_LB),
      CI_UB  = cross(bounds$CI_UB)
    )
  )
}

## ---- Sensitivity plot function ------------------------------
#' Sensitivity Plot
#'
#' Plot ATT bounds and confidence intervals as a function of the sensitivity
#' parameter \eqn{\log \Gamma} (Park & Tchetgen Tchetgen, 2026+).
#'
#' @param output Output from \code{\link{UDID_Nonparametric}},
#'   \code{\link{UDID_Parametric}}, or \code{\link{UDID_Sensitivity_Bounds}}.
#' @param alpha Confidence level (default 0.95).
#' @param col_att_fill Fill color for ATT band.
#' @param col_att_line Line color for ATT bounds.
#' @param col_ci_fill Fill color for CI band.
#' @param col_ci_line Line color for CI bounds.
#' @param xlab X-axis label.
#' @param ylab Y-axis label.
#' @param main Plot title.
#' @param zero_line Logical; draw horizontal line at zero.
#' @param legend Logical; show legend.
#' @param legend.pos Legend position.
#'
#' @details
#' \code{UDID_Sensitivity_Plot} visualizes the sensitivity analysis results
#' by plotting:
#' \itemize{
#'   \item A shaded dark band for the ATT upper and lower bounds
#'   (\code{ATT_UB}, \code{ATT_LB}) as a function of \eqn{\log \Gamma}.
#'   \item A shaded light band for the pointwise confidence interval bounds
#'   (\code{CI_UB}, \code{CI_LB}).
#'   \item Vertical dashed lines at zero-crossing points, indicating the
#'   values of \eqn{\log \Gamma} at which the ATT bounds or CI bounds cross
#'   zero.
#' }
#' A horizontal reference line at zero is drawn by default.
#' The zero-crossing values indicate how much unmeasured confounding
#' (measured in \eqn{\log \Gamma}) would be required to either nullify
#' the treatment effect estimate or render it statistically insignificant.
#'
#' @return Invisibly returns the crossings list from
#'   \code{\link{UDID_Sensitivity_Bounds}}.
#'
#' @references
#' \itemize{
#'   \item Park, C., & Tchetgen Tchetgen, E. (2026+).
#'     A Universal Nonparametric Framework for Difference-in-Differences Analyses.
#'     \url{https://arxiv.org/abs/2212.13641}.
#' }
#'
#' @seealso \code{\link{UDID_Sensitivity_Bounds}} for the underlying
#'   computation, \code{\link{UDID_Nonparametric}} and
#'   \code{\link{UDID_Parametric}} for the estimators.
#'
#' @export
UDID_Sensitivity_Plot <- function(output,
                                  alpha        = 0.95,
                                  col_att_fill = grDevices::rgb(0, 0, 0, 0.4),
                                  col_att_line = grDevices::rgb(0, 0, 0, 0.75),
                                  col_ci_fill  = grDevices::rgb(1, 0, 0, 0.25),
                                  col_ci_line  = grDevices::rgb(1, 0, 0, 0.6),
                                  xlab         = expression("log " * Gamma),
                                  ylab         = "ATT",
                                  main         = "Sensitivity Analysis",
                                  zero_line    = TRUE,
                                  legend       = TRUE,
                                  legend.pos   = "bottomleft") {
  
  sb <- UDID_Sensitivity_Bounds(output, alpha)
  result <- sb$bounds
  
  x     <- result$log_Gamma
  ylim  <- range(c(result$CI_LB, result$CI_UB), na.rm = TRUE)
  ylim  <- ylim + c(-1, 1) * diff(ylim) * 0.05
  
  graphics::plot(x, result$ATT_LB, type = "n", ylim = ylim,
       xlab = xlab, ylab = ylab, main = main)
  
  ## CI band
  graphics::polygon(c(x, rev(x)), c(result$CI_LB, rev(result$CI_UB)),
          col = col_ci_fill, border = NA)
  graphics::lines(x, result$CI_LB, col = col_ci_line)
  graphics::lines(x, result$CI_UB, col = col_ci_line)
  
  ## ATT band
  graphics::polygon(c(x, rev(x)), c(result$ATT_LB, rev(result$ATT_UB)),
          col = col_att_fill, border = NA)
  graphics::lines(x, result$ATT_LB, col = col_att_line)
  graphics::lines(x, result$ATT_UB, col = col_att_line)
  
  if (zero_line)
    graphics::abline(h = 0, lty = 3, col = "grey40")
  
  label_cross <- function(v, col) {
    if (!is.na(v)) {
      graphics::abline(v = v, lty = 2, col = col)
      graphics::text(x = v, y = ylim[1], labels = round(v, 3),
           col = col, adj = c(0.5, -0.4), cex = 0.8, pos=2)
    }
  }
  label_cross(sb$crossings$ATT_LB, col_att_line)
  label_cross(sb$crossings$ATT_UB, col_att_line)
  label_cross(sb$crossings$CI_LB,  col_ci_line)
  label_cross(sb$crossings$CI_UB,  col_ci_line)
  
  if (legend) {
    graphics::legend(legend.pos, bty = "n",
           legend = c("ATT bounds",
                      paste0(round(alpha * 100), "% CI bounds")),
           fill   = c(col_att_fill, col_ci_fill),
           border = c(col_att_line, col_ci_line))
  }
  
  invisible(sb$crossings)
}
