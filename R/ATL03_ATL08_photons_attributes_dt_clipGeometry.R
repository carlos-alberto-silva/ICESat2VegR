#' Clip Joined ATL03 and ATL08 by Geometry
#'
#' @description This function clips joined ATL03 and ATL08 photon attributes within a given geometry.
#'
#' @param atl03_atl08_dt An S4 object of class [`ICESat2VegR::icesat2.atl03atl08_dt-class`] containing ATL03 and ATL08 data
#' (output of [`ICESat2VegR::ATL03_ATL08_photons_attributes_dt_join()`] function).
#' @param polygon An object of class [`terra::SpatVector`], which can be loaded as an ESRI shapefile using [terra::vect] function.
#' @param split_by Optional. Polygon ID. If defined, the data will be clipped by each polygon using
#' the polygon ID from the attribute table.
#'
#' @return Returns an S4 object of class [`ICESat2VegR::icesat2.atl03atl08_dt-class`]
#' containing the clipped ATL08 attributes.
#'
#' @examples
#' # Specifying the path to ATL03 and ATL08 files
#' atl03_path <- system.file("extdata", "atl03_clip.h5", package = "ICESat2VegR")
#' atl08_path <- system.file("extdata", "atl08_clip.h5", package = "ICESat2VegR")
#'
#' # Reading ATL03 and ATL08 data (h5 files)
#' atl03_h5 <- ATL03_read(atl03_path)
#' atl08_h5 <- ATL08_read(atl08_path)
#'
#' # Joining ATL03 and ATL08 photon attributes
#' atl03_atl08_dt <- ATL03_ATL08_photons_attributes_dt_join(atl03_h5, atl08_h5)
#' head(atl03_atl08_dt)
#'
#' # Specifying the path to the shapefile
#' polygon_filepath <- system.file("extdata", "clip_geom.shp", package = "ICESat2VegR")
#'
#' # Reading shapefile as a SpatVector object
#' polygon <- terra::vect(polygon_filepath)
#'
#' # Clipping ATL08 terrain attributes by geometry
#' atl03_atl08_dt_clip <- ATL03_ATL08_photons_attributes_dt_clipGeometry(
#'   atl03_atl08_dt,
#'   polygon,
#'   split_by = "id"
#' )
#' head(atl03_atl08_dt_clip)
#'
#' close(atl03_h5)
#' close(atl08_h5)
#' @export
ATL03_ATL08_photons_attributes_dt_clipGeometry <- function(atl03_atl08_dt, polygon, split_by = NULL) {
  if (!inherits(atl03_atl08_dt, "icesat2.atl03atl08_dt")) {
    stop("atl03_atl08_dt must be an object of class 'icesat2.at03atl08_dt'")
  }

  nid <- NA

  exshp <- terra::ext(polygon)

  # Clip data by bounding box of the polygon
  atl03_atl08_dt <- ATL03_ATL08_photons_attributes_dt_clipBox(
    atl03_atl08_dt,
    lower_left_lon = exshp$xmin,
    upper_right_lon = exshp$xmax,
    lower_left_lat = exshp$ymin,
    upper_right_lat = exshp$ymax
  )

  # Remove NA values
  atl03_atl08_dt <- na.omit(atl03_atl08_dt)

  if (nrow(atl03_atl08_dt) == 0) {
    warning("The polygon does not overlap the ATL08 data")
    empty <- data.table()
    prepend_class(empty, "icesat2.atl03atl08_dt")
    return(empty)
  }

  atl03_atl08_dt[, nid := .I]
  points <- to_vect(atl03_atl08_dt)
  points$rowNumber <- seq_along(points)
  pts <- terra::intersect(terra::makeValid(points), terra::makeValid(polygon))

  if (!is.null(split_by)) {
    if (split_by %in% names(polygon)) {
      newFile <- atl03_atl08_dt[pts$nid, ]
      newFile$poly_id <- pts[[split_by]]
    } else {
      stop("The ", split_by, " is not included in the attribute table. Please check the names in the attribute table.")
    }
  } else {
    newFile <- atl03_atl08_dt[pts$nid, ]
  }

  prepend_class(newFile, "icesat2.atl03atl08_dt")
  return(newFile)
}
