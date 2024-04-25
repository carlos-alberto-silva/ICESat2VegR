#' @include model_tools.R

setRefClass("ee.Classifier")

#' `predict` method for generating predictions for the ee.Classifier.
#'
#' @param x The ee.Classifier generated either by using the regular
#' `ee` API functions or using the package's [`build_ee_forest()`] function.
#' @param data [`data.frame-class`] or `ee.FeatureCollection` with the
#' features input to be used for the prediction.
#'
#' @return A [`numeric-class`] with the predicted values
#'
#' @export
"predict.ee.Classifier" <- function(x, data) {
  if (!inherits(data, "ee.featurecollection.FeatureCollection")) {
    data <- build_fc(data, NULL)
  }

  predicted <- data$classify(x)
  classification <- predicted$aggregate_array("classification")$getInfo()
  return(classification)
}


"predict.yai" <- function(object, ...) {
  x <- list(...)[[1]]
  rownames(x) = paste0("_", rownames(x))
  y <- object$yRefs
  new.t <- yaImpute::newtargets(object, newdata = cbind(x))
  
  res <- yaImpute::impute(new.t,vars=yaImpute::yvars(new.t))
  rownames(res) <- NULL
  res[,!apply(res, 2, anyNA)]
}
