#' Sensitivity Bounds
#'
#' Compute ATT and CI bounds from a sensitivity analysis output.
#'
#' @param output A data.frame with columns \code{log_Gamma}, \code{ATT}, \code{SE}.
#' @param alpha Confidence level (default 0.95).
#'
#' @return A list with \code{bounds} (data.frame) and \code{crossings} (list of zero-crossing values).
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
#' Plot ATT and confidence interval as a function of the sensitivity parameter.
#'
#' @param output A data.frame (or list with a \code{bounds} element) with
#'   columns \code{log_Gamma}, \code{ATT}, \code{SE}.
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
#' @return Invisibly returns the crossings list.
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
