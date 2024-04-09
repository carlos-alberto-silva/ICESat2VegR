#' Applies the model in the input rasters to generate the
#' modelled raster
#'
#' @param model [`randomForest::randomForest-class`]. The classifier (regression) model
#' to use for the prediction.
#' @param stack. A vector/list of ee.Image created using the ee API.
#'
#' @return An ee.Image resulting from the model.
#'
#' @export
map_create <- function(model, stack) {
  trees <- build_forest(model)
  gee_rf <- ee$Classifier$decisionTreeEnsemble(trees)
  result <- stack$classify(gee_rf)
  result
}
