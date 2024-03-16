#' Based on observed and predicted values create analyzes the model
#' accuracy based on multiple statistics
#' 
#' @param obs [`numeric-class`]. A numeric vector for the observed values.
#' @param est [`numeric-class`]. A numeric vector for predicted values.
#' @param plotR [`logical-class`]. If TRUE will plot the observed x predicted
#' 1 to 1 plot.
#' @param xlab [`character-class`]. The label for x axis, default "Observed".
#' @param ylab [`character-class`]. The label for y axis, default "Predicted".
#' 
#' @export
stat_model <- function(obs, est , plotR=TRUE, xlab="Observed", ylab="Predicted", ...)
{

  xy<-na.omit(cbind(obs,est))

  obs<-xy[,1]
  est<-xy[,2]

  library(formattable)

  rmse <- sqrt( sum(( est - obs )^2)/length(obs) ) # Root mean square error
  rmse <-formattable(rmse, digits = 2, format = "f")
  bias <- mean( est - obs ) # bias
  bias<-formattable(bias, digits = 2, format = "f")
  rmseR <- 100 * sqrt( sum(( est - obs )^2)/length(obs) ) / mean( obs )
  biasR <- 100 * mean( est - obs ) / mean( obs )
  MAE<- sum((  est - obs)^2)/length(obs) # Mean absolute error
  PMRE<- sum(sqrt((( est - obs )/obs)^2)/length(obs))*100 # percent mean relative error. (PMRE)


  r <- cor(est,obs)
  lms<-lm(obs~est)
  r2<-summary(lms)$r.squared
  StatInfo<-data.frame( Stat=c("rmse","rmseR","bias","biasR","MAE","PMRE","r","r2","intercept","slope"),
                        Values=round(c(rmse,rmseR,bias,biasR,MAE,PMRE,r,r2,coef(lms)[1],coef(lms)[2]),2)) #0 valor 6 eh o numero de casas decimais depois do ponto.


  StatInfo2<-StatInfo
  StatInfo2$Values<-round(StatInfo2$Values,2)
  if (plotR==TRUE) {
    plot(obs, est, xlab=xlab, ylab=ylab, ...)
    abline(0,1, col="red")
    abline(lm(est~obs), col="black", lwd=2)
    legend("topleft",c(paste0("RMSE=",StatInfo2[1,2],"Mg/ha"),
                       paste0("RMSE=",StatInfo2[2,2]," %"),
                       paste0("Bias=",StatInfo2[3,2],"Mg/ha"),
                       paste0("Bias=",StatInfo2[4,2]," %"),
                       paste0("adj.R2=",StatInfo2[8,2]),
                       paste0("Inter=",StatInfo2[9,2]),
                       paste0("Slope=",StatInfo2[10,2])),bty="n")
  }

  return(StatInfo)
}