#' @export
setGeneric(
  "addEEImage",
  def = function(map, x, bands, min_value = 0, max_value = 1, palette = c("red", "green"), ...) {
    standardGeneric("addEEImage")
  }
)


setRefClass("leaflet")
#' Adds Earth Engine Image class to leaflet
#'
#' @param map [`leaflet-class`]. A leaflet class map.
#' @param x Earth Engine Image open with `ee$Image`.
#' @param min_value The minimum value to represent for visualization purposes
#' @param max_value The maximum value to represent for visualization purposes
#' @param palette A [`character-class`] describing a list of colors in either
#' hexadecimal or valid CSS color names.
#' @param ... Other parameters to be passed to [`leaflet::addTiles()`] function.
#'
#' @return A [`leaflet-class`]/[`htmlwidget-class`] with the Earth Engine image added.
#'
#' @examples
#' collection_id <- "NASA/HLS/HLSL30/v002"
#'
#' ee_collection <- ee$ImageCollection(collection_id)
#'
#' cloudMask <-
#'   2^1 +
#'   2^2 +
#'   2^3
#'
#' hlsMask <- function(image) {
#'   image$updateMask(
#'     !(image[["Fmask"]] & cloudMask)
#'   )
#' }
#'
#' image <- c(
#'   ee_collection$filterDate("2020-04-01", "2020-07-31")$map(hlsMask),
#'   ee_collection$filterDate("2021-04-01", "2021-07-31")$map(hlsMask)
#' )$filter("CLOUD_COVERAGE < 30")$median()
#'
#'
#' library(leaflet)
#' leaflet() %>%
#'   addEEImage(
#'     image,
#'     bands = list(3, 2, 1),
#'     min_value = 0.001,
#'     max_value = 0.2
#'   ) %>%
#'   setView(lng = -82.2345, lat = 29.6552, zoom = 10)
#' @include gee-base.R
#' @export
setMethod(
  "addEEImage",
  c("leaflet", "ee.image.Image"),
  function(map, x, bands, min_value = 0, max_value = 1, palette = c("#d55e00", "#cc79a7", "#f0e442", "#0072b2", "#009e73"), ...) {
    stopifnot("`bands` list should be either one or three bands" = (length(bands) == 1 || length(bands) == 3))

    if (length(bands) == 3) {
      vis <- x[[bands]]$visualize(
        min = min_value,
        max = max_value
      )
    } else {
      vis <- x[[bands]]$visualize(
        min = min_value,
        max = max_value,
        palette = palette
      )
    }


    url <- vis$getMapId()$tile_fetcher$url_format
    print(url)

    return(
      leaflet::addTiles(
        map,
        urlTemplate = url,
        ...
      )
    )
  }
)
