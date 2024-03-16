#' @export
map_create <- function(model, stack) {
  result <- stack$classify(model)
  result
}