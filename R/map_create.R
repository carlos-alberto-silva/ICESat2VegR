#' Applies the model in the input rasters to generate the
#' modelled raster
#' 
#' @param model ee.Classfier. The classifier (regression) model
#' to use for the prediction.
#' @param stack. A vector/list of ee.Image created using the ee API.
#' 
#' @return An ee.Image resulting from the model.
#' 
#' @export
map_create <- function(model, stack) {
  result <- stack$classify(model)
  result
}