#' Generic function to export icesat2 classes to [`terra::SpatVector-class`]
#' @param x The icesat2 object to convert to [`terra::SpatVector-class`]
#' @param ... other potential parameters needed.
#'
#' @return The icesat2 object as the appropriate [`terra::SpatVector-class`]
#'
#' @examples
#' # Get path to ATL03 h5 file
#' atl03_path <- system.file("extdata",
#'   "atl03_clip.h5",
#'   package = "ICESat2VegR"
#' )
#' # Reading ATL03 data (h5 file)
#' atl03_h5 <- ATL03_read(atl03_path = atl03_path)
#' # Extracting ATL03 segment attributes
#' atl03_segment_dt <- ATL03_seg_attributes_dt(atl03_h5 = atl03_h5)
#'
#' # Extract vector
#' atl03_segment_vect <- to_vect(atl03_segment_dt)
#'
#' # Terra plot
#' terra::plot(atl03_segment_vect, col = atl03_segment_vect$segment_ph_cnt)
#'
#' \donttest{
#' # Export as temp gpkg
#' temp_vect <- tempfile(fileext = ".gpkg")
#' terra::writeVector(atl03_segment_vect, temp_vect)
#' }
#'
#' head(atl03_segment_dt)
#' close(ATL03_h5)
#' @include class.icesat2.R
#' @export
setGeneric(
  "to_vect",
  def = function(x, ...) {
    standardGeneric("to_vect")
  }
)

generic_to_vect <- function(x, lon, lat) {
  df <- as.data.frame(x)
  terra::vect(
    df,
    geom = c(lon, lat),
    crs = "epsg:4326"
  )
}


#' @export
setMethod(
  "to_vect",
  "icesat2.atl03_seg_dt",
  function(x, ...) {
    generic_to_vect(x, "reference_photon_lon", "reference_photon_lat")
  }
)


#' @export
setMethod(
  "to_vect",
  "icesat2.atl08_dt",
  function(x, ...) {
    generic_to_vect(x, "longitude", "latitude")
  }
)

#' @export
setMethod(
  "to_vect",
  "icesat2.atl03atl08_dt",
  function(x, ...) {
    generic_to_vect(x, "lon_ph", "lat_ph")
  }
)

#' @export
setMethod(
  "to_vect",
  "icesat2.atl03_dt",
  function(x, ...) {
    generic_to_vect(x, "lon_ph", "lat_ph")
  }
)

#' @export
setMethod(
  "to_vect",
  "icesat2.atl03_atl08_seg_dt",
  function(x, ...) {
    generic_to_vect(x, "longitude", "latitude")
  }
)
