#' Add an Earth Engine Image to a leaflet map
#'
#' @param map A `leaflet::leaflet` widget.
#' @param img An Earth Engine image (reticulate `python.builtin.object` with `$visualize`, `$select`).
#' @param bands Character or numeric vector (length 1 or 3). If `NULL`, uses all bands from the image.
#' @param min_value Numeric. Minimum visualization range.
#' @param max_value  Numeric. Maximum visualization range.
#' @param palette Character vector of colors for single-band rendering.
#' @param group Optional overlay group name.
#' @param aoi Optional EE geometry to clip before rendering.
#' @param ... Passed to `leaflet::addTiles()` (e.g., `options = leaflet::tileOptions(opacity = 0.8)`).
#'
#' @return The input `leaflet` map with the EE image layer added.
#'
#' @examples
#' \dontrun{
#' # Initialize Earth Engine (use your own project id)
#' Sys.setenv(EE_PROJECT = "your-ee-project-id")
#' ee_initialize()
#'
#' # Harmonized Landsat-Sentinel (HLS) surface reflectance collection
#' collection_id <- "NASA/HLS/HLSL30/v002"
#' ee_collection <- ee$ImageCollection(collection_id)
#'
#' # Mask cloud, cloud-shadow and adjacent-cloud pixels (Fmask bits 1-3)
#' cloudMask <- 2^1 + 2^2 + 2^3
#' hlsMask <- function(image) {
#'   image$updateMask(!(image[["Fmask"]] & cloudMask))
#' }
#'
#' # Cloud-free median composite over two summer seasons
#' image <- c(
#'   ee_collection$filterDate("2020-04-01", "2020-07-31")$map(hlsMask),
#'   ee_collection$filterDate("2021-04-01", "2021-07-31")$map(hlsMask)
#' )$filter("CLOUD_COVERAGE < 30")$median()
#'
#' # Show the RGB composite (bands 3-2-1) on an interactive map
#' if (require("leaflet")) {
#'   leaflet() %>%
#'     addEEImage(
#'       image,
#'       bands = c(3, 2, 1),
#'       min_value = 0.001,
#'       max_value = 0.2
#'     ) %>%
#'     setView(lng = -82.2345, lat = 29.6552, zoom = 10)
#' }
#' }
addEEImage <- function(map,
                       img,
                       bands = NULL,
                       min_value = 0,
                       max_value = 1,
                       palette = c("#d55e00", "#cc79a7", "#f0e442", "#0072b2", "#009e73"),
                       group = NULL,
                       aoi = NULL,
                       ...) {
  .must_have("leaflet"); .must_have("reticulate")

  # allow min/max aliases via ..., but don't use missing() on list elements
  dots <- list(...)
  if ("min" %in% names(dots)) min_value <- dots[["min"]]
  if ("max" %in% names(dots)) max_value <- dots[["max"]]

  .add_ee_image_layer(
    map       = map,
    x         = img,
    bands     = bands,
    aoi       = aoi,
    group     = group,
    min_value = min_value,
    max_value = max_value,
    palette   = palette,
    is_class  = FALSE,
    scale_to_int = TRUE,
    int_factor   = 1000L
  )
}
