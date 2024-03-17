#' Given a stack image raster from GEE
#' retrieve the point geometry with values for the images
#'
#' @param stack A single image or a vector/list of images from Earth Engine.
#' @param geom A geometry from [`terra::SpatVector-class`] read with [`terra::vect`].
#' @param scale The scale in meters for the extraction (image resolution).
#' 
#' @return A [data.table::data.table-class] with the properties extracted
#' from the ee images.
#'
#' @export
seg_gee_ancillary_dt_extract <- function(stack, geom, scale = 30, chunk_size = 10) {
  n <- nrow(geom)
  dts <- list()
  for (ii in seq(1, n, chunk_size)) {
    tail <- ii + chunk_size - 1
    if (tail > n) {
      tail <- n
    }
    message(gettextf("Processing %d-%d of %d", ii, tail, n))
    sampled <- extract(stack, geom[ii:tail], scale)
    dts[[""]] <- ee_to_dt(sampled)
  }
  dt <- data.table::rbindlist(dts)
  return (dt)
}


#' Given a geometry with point samples and images from Earth Engine
#' retrieve the point geometry with values for the images
#'
#' @param stack A single image or a vector/list of images from Earth Engine.
#' @param geom A geometry from [`terra::SpatVector-class`] read with [`terra::vect()`].
#' @param scale The scale in meters for the extraction (image resolution).
#' 
#' @return An ee.FeatureCollection with the properties extracted from the 
#' stack of images from ee. 
#'
#' @export
extract <- function(stack, geom, scale) {
  final <- c(stack)
  tempjson <- tempfile(fileext = ".geojson")
  terra::writeVector(geom, tempjson, filetype = "geojson")
  parsed <- jsonlite::parse_json(readLines(tempjson))
  geojson <- ee$FeatureCollection(parsed)

  sampled <- final$sampleRegions(
    collection = geojson,
    scale = scale,
    tileScale = 16
  )

  return(sampled)
}