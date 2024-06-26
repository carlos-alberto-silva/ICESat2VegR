#' Adds Earth Engine Image class to leaflet
#'
#' @param map `leaflet::leaflet`. A leaflet class map.
#' @param x Earth Engine Image open with `ee$Image`.
#' @param bands define the bands that should be used for visualization, default NULL,
#' will use the available bands from the Image.
#' @param min_value The minimum value to represent for visualization purposes
#' @param max_value The maximum value to represent for visualization purposes
#' @param palette A [`character-class`] describing a list of colors in either
#' hexadecimal or valid CSS color names.
#' @param ... Other parameters to be passed to `leaflet::addTiles()` function.
#'
#' @return A `leaflet::leaflet`/`htmlwidget` with the Earth Engine image added.
#'
#' @examples
#' \dontrun{
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
#' if (require("leaflet")) {
#'   leaflet() %>%
#'     addEEImage(
#'       image,
#'       bands = list(3, 2, 1),
#'       min_value = 0.001,
#'       max_value = 0.2
#'     ) %>%
#'     setView(lng = -82.2345, lat = 29.6552, zoom = 10)
#'   }
#' }
#' @include gee-base.R
#' @export
setGeneric(
  "addEEImage",
  def = function(map, x, bands, min_value = 0, max_value = 1, palette = c("red", "green"), ...) {
    if (requireNamespace("leaflet", quietly = TRUE)) {
      standardGeneric("addEEImage")
    } else {
      stop("The 'leaflet' package is required to use this function.")
    }
  }
)

defaultPallete <- c("#d55e00", "#cc79a7", "#f0e442", "#0072b2", "#009e73")
setRefClass("leaflet")

#' Adds Earth Engine Image class to leaflet map
#'
#' This method handles adding an Earth Engine Image to a leaflet map.
#'
#' @param map A `leaflet::leaflet` object. A leaflet map to which the Earth Engine image will be added.
#' @param x An Earth Engine Image object created with `ee$Image`.
#' @param bands A vector defining the bands to be used for visualization. Default is NULL, which uses the available bands from the image.
#' @param min_value The minimum value for visualization purposes. Default is 0.
#' @param max_value The maximum value for visualization purposes. Default is 1.
#' @param palette A character vector describing a list of colors in either hexadecimal or valid CSS color names. Default is a predefined palette.
#' @param ... Other parameters to be passed to the `leaflet::addTiles()` function.
#'
#' @return A `leaflet::leaflet`/`htmlwidget` object with the Earth Engine image added.
#' @export
setMethod(
  "addEEImage",
  c("leaflet", "ee.image.Image"),
  function(map, x, bands = NULL, min_value = 0, max_value = 1, palette = defaultPallete, ...) {
    if (is.null(bands)) {
      bands <- names(x)
    }

    stopifnot(
      "Image should have either one or three bands, if not define the `bands` parameter with either 1 or 3 bands" =
        (length(bands) == 1 || length(bands) == 3)
    )

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

    return(
      invisible(
        leaflet::addTiles(
          map,
          urlTemplate = url,
          ...
        )
      )
    )
  }
)
