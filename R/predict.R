#' @include model_tools.R

setRefClass("ee.Classifier")

#' `predict` method for generating predictions for the ee.Classifier.
#' 
#' @param x. The ee.Classifier generated either by using the regular
#' `ee` API functions or using the package's [`model_fit()`] function.
#' @param data [`data.frame-class`] or `ee.FeatureCollection` with the
#' features input to be used for the prediction.
#' 
#' @return A [`numeric-class`] with the predicted values
#' 
#' @export
"predict.ee.Classifier" <- function(x, data, ...) {
  if (!inherits(data, "ee.featurecollection.FeatureCollection")) {
    data <- build_fc(data, NULL)
  }
  
  predicted <- data$classify(x)
  classification <- predicted$aggregate_array("classification")$getInfo()
  return(classification)
}
