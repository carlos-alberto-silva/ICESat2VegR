#' Map create
#'
#' @description Applies the fitted model as [`randomForest::randomForest`] object to
#' the input raster stack (as ee.Image) to produce a wall-to-wall maps.
#'
#' @param model [`randomForest::randomForest`]. The classifier (regression) model
#' to use for the prediction.
#' @param stack A vector/list of ee.Image created using the ee API.
#'
#' @return An `ee.Image` resulting from applying the fitted model.
#'
#' @export
map_create <- function(model, stack) {
  trees <- build_ee_forest(model)
  gee_rf <- ee$Classifier$decisionTreeEnsemble(trees)
  result <- stack$classify(gee_rf)
  result
}
