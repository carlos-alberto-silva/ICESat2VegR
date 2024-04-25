#' Model fit stats
#'
#' @description Computes absolute and relative root-mean-square error (RMSE) and bias, and adjusted coefficient of determination
#' from a linear relationship between predicted and observed data
#'
#' @param observed  A vector containing observed data
#' @param predicted A vector containing predicted data
#' @param plotstat if TRUE, a plot showing the 1:1 figure will be displayed. Default is TRUE
#' @param legend Character for legend position. Default is "topleft"
#' @param unit Character indicating the unit for the observed and predicted data (e.g."Mg/ha")
#' @param ... Other params to be redirected to the plot.
#'
#' @return Returns an data.frame object containing the list of stats
#'
#' @examples
#'
#'Observed and predicted aboveground biomass
#'observed<-c(178,33,60,80,104,204,146)
#'predicted<-c(184,28.5,55,85,105,210,155)
#'
#'stats_model(observed, predicted , plotstat=TRUE, legend="topleft", unit="Mg/ha",
#'            xlab="Obserbed AGBD (Mg/ha)",
#'            ylab="Predicted AGBD (Mg/ha)", pch=16)
#'
#' @export
stats_model <- function(observed, predicted, plotstat = TRUE, legend = "topleft", unit = " m", ...) {
  rmse <- sqrt(sum((predicted - observed)^2) / length(observed)) # Root mean square error
  bias <- mean(predicted - observed) # bias
  rmseR <- 100 * sqrt(sum((predicted - observed)^2) / length(observed)) / mean(observed)
  biasR <- 100 * mean(predicted - observed) / mean(observed)

  r <- cor(predicted, observed)
  r2 <- summary(lm(observed ~ predicted))$r.squared
  StatInfo <- data.frame(
    Stat = c("rmse", "rmseR", "bias", "biasR", "r", "r2"),
    Values = round(c(rmse, rmseR, bias, biasR, r, r2), 10)
  )

  # StatInfo2<-StatInfo
  # StatInfo2$Values<-round(StatInfo2$Values,2)
  if (plotstat == TRUE) {
    plot(observed, predicted, ...)
    abline(0, 1, col = "red", lwd = 2)
    abline(lm(predicted ~ observed), col = "black", lwd = 2)
    legend("topleft", c(
      gettextf("RMSE=%.4f %s", StatInfo[1, 2], unit),
      gettextf("RMSE=%.4f%s", StatInfo[2, 2], " %"),
      gettextf("Bias=%.4f%s", StatInfo[3, 2], unit),
      gettextf("Bias=%.4f%s", StatInfo[4, 2], " %"),
      gettextf("adj.R2=%.4f", StatInfo[6, 2])
    ), bty = "n")
  }

  return(StatInfo)
}
