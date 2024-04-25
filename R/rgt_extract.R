#' Extract reference ground track from ATL03 segments
#'
#' @description This function extracts reference ground track from ICESat-2 ATL03 data
#' and writes it as a GDAL vector format.
#'
#' @param h5 A ICESat-2 ATL03 object (output of [`ATL03_read()`] or [`ATL08_read()`] function).
#' An S4 object of class [`ICESat2VegR::icesat2.atl03_dt-class`] or [`ICESat2VegR::icesat2.atl08_dt-class`].
#'
#' @return Returns the ground track boundaries as [`terra::SpatVector-class`] extracted from "orbit_info",
#' along with other orbit information.
#'
#' @details This function will use the reference photons from the segments as reference
#' for deriving the ground tracks. The begining and end of the lines are interpolated from
#' the information regarding the position of the reference photon within the segment and the
#' segment length.
#'
#' @seealso \url{https://icesat-2.gsfc.nasa.gov/sites/default/files/page_files/ICESat2_ATL03_ATBD_r006.pdf}
#'
#' @examples
#'
#' # Specifying the path to ATL03 H5 file
#' atl03_path <- system.file("extdata",
#'   "atl03_clip.h5",
#'   package = "ICESat2VegR"
#' )
#'
#' # Reading ATL03 data (h5 file)
#' atl03_h5 <- ATL03_read(atl03_path = atl03_path)
#'
#' # Extracting ATL03 photons attributes
#' rgt <- rgt_extract(h5 = atl03_h5)
#' head(rgt)
#' 
#' terra::plet(rgt)
#'
#' close(atl03_h5)
#' @include class.icesat2.R
#' @importFrom data.table data.table rbindlist
#' @importFrom terra writeVector
#' @export
rgt_extract <- function(h5) {
  stopifnot(
    "The input file needs to be one ATL03_read or ATL08_read products" = inherits(h5, "icesat2.h5")
  )

  rgt_data <- data.frame(
    cycle_number = h5[["orbit_info/cycle_number"]][],
    orbit_number = h5[["orbit_info/orbit_number"]][],
    rgt = h5[["orbit_info/rgt"]][]
  )

  rgt <- terra::vect(
    cbind(
      h5[["orbit_info/bounding_polygon_lon1"]][],
      h5[["orbit_info/bounding_polygon_lat1"]][]),
    type = "polygons",
    atts = rgt_data,
    crs = "epsg:4326"
  )

  return(rgt)
}
