#' Model fit metrics
#'
#' @description Computes RMSE (absolute and relative), MAE (absolute and relative),
#' bias (absolute and relative), Pearson correlation (r), and adjusted R² from a
#' linear fit between predicted and observed values. Optionally draws a 1:1 plot
#' with the regression line.
#'
#' @param observed  Numeric vector of observed values.
#' @param predicted Numeric vector of predicted values (same length/order as `observed`).
#' @param plotstat  Logical; if TRUE, draws a 1:1 plot with the regression line. Default: FALSE.
#' @param legend    Character position for the plot legend (e.g., "topleft"). Default: "topleft".
#' @param unit      Character indicating the unit for RMSE/MAE/Bias (e.g., "Mg/ha"). Default: "".
#' @param ...       Additional arguments passed to [graphics::plot()].
#'
#' @return A data.frame with rows for `rmse`, `rmseR` (%), `mae`, `maeR` (%),
#' `bias`, `biasR` (%), `r`, and `adj_r2`.
#' @examples
#' observed  <- c(178, 33, 60, 80, 104, 204, 146)
#' predicted <- c(184, 28.5, 55, 85, 105, 210, 155)
#' fit_metrics(observed, predicted,
#'             plotstat = TRUE, legend = "topleft", unit = "Mg/ha",
#'             xlab = "Observed AGBD (Mg/ha)", ylab = "Predicted AGBD (Mg/ha)", pch = 16)
#' @export
fit_metrics <- function(observed, predicted,
                        plotstat = FALSE, legend = "topleft", unit = "", ...) {
  # coerce and align
  observed  <- as.numeric(observed)
  predicted <- as.numeric(predicted)
  if (length(observed) != length(predicted)) {
    stop("`observed` and `predicted` must have the same length.")
  }
  ok <- stats::complete.cases(observed, predicted)
  if (!all(ok)) {
    observed  <- observed[ok]
    predicted <- predicted[ok]
  }
  n <- length(observed)
  if (n < 2L) stop("Need at least 2 complete pairs to compute statistics.")

  # core errors
  err   <- predicted - observed
  rmse  <- sqrt(mean(err^2))
  mae   <- mean(abs(err))
  bias  <- mean(err)

  # relative (%) metrics w/ safe denominator
  m_obs <- mean(observed)
  if (isTRUE(all.equal(m_obs, 0))) {
    rmseR <- NA_real_
    maeR  <- NA_real_
    biasR <- NA_real_
  } else {
    rmseR <- 100 * rmse / m_obs
    maeR  <- 100 * mae  / m_obs
    biasR <- 100 * bias / m_obs
  }

  # association and model fit
  r <- suppressWarnings(stats::cor(predicted, observed))
  if (!is.finite(r)) r <- NA_real_

  fit    <- stats::lm(observed ~ predicted)
  adj_r2 <- summary(fit)$adj.r.squared

  out <- data.frame(
    stat  = c("rmse", "rmseR", "mae", "maeR", "bias", "biasR", "r", "adj_r2"),
    value = c(rmse,   rmseR,   mae,   maeR,   bias,   biasR,   r,  adj_r2),
    unit  = c(unit,   "%",     unit,  "%",    unit,   "%",    "", "")
  )

  if (isTRUE(plotstat)) {
    graphics::plot(observed, predicted, ...)
    graphics::abline(0, 1, col = "red", lwd = 2)   # 1:1 line
    graphics::abline(fit, col = "black", lwd = 2)  # regression line
    lab <- c(
      sprintf("RMSE = %.4f %s", out$value[out$stat=="rmse"],  unit),
      sprintf("RMSE = %.4f%%",   out$value[out$stat=="rmseR"]),
      sprintf("MAE  = %.4f %s", out$value[out$stat=="mae"],   unit),
      sprintf("MAE  = %.4f%%",   out$value[out$stat=="maeR"]),
      sprintf("Bias = %.4f %s", out$value[out$stat=="bias"],  unit),
      sprintf("Bias = %.4f%%",   out$value[out$stat=="biasR"]),
      sprintf("adj.R² = %.4f",   out$value[out$stat=="adj_r2"])
    )
    graphics::legend(legend, lab, bty = "n")
  }

  out
}

