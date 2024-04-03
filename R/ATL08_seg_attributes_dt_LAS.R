#' Converts ATL08 segments to LAS
#'
#' @param atl08_dt  An S4 object of class [ICESat2VegR::icesat2.atl08_dt] containing ATL03 and ATL08 data
#' (output of [ATL03_ATL08_photons_attributes_dt_join()] function).
#' @param output character. The output path of for the LAS(Z) file(s)
#' The function will create one LAS file per UTM Zone in WGS84 datum.
#'
#' @return Nothing, it just saves outputs as LAS file in disk
#'
#' @details
#' As the las format expects a metric coordinate reference system (CRS)
#' we use helper functions to define UTM zones to which the original
#' ICESat-2 data will be converted.
#'
#' The function credits go to Chuck Gantz- chuck.gantz@globalstar.com.
#'
#' @seealso
#' https://oceancolor.gsfc.nasa.gov/docs/ocssw/LatLong-UTMconversion_8cpp_source.html
#'
#' @examples
#' # Specifying the path to ATL08 file
#' atl08_path <- system.file("extdata",
#'   "atl08_clip.h5",
#'   package = "ICESat2VegR"
#' )
#'
#' # Reading ATL08 data (h5 file)
#' atl08_h5 <- ATL08_read(atl08_path = atl08_path)
#'
#' # # Extracting ATL03 and ATL08 photons and heights
#' atl08_dt <- ATL08_seg_attributes_dt(atl08_h5)
#'
#' outputLaz <- tempfile(fileext = ".laz")
#' ATL08_seg_attributes_dt_LAS(
#'   atl08_dt,
#'   outputLaz
#' )
#'
#' close(atl08_h5)
#' @include lasTools.R
#' @importFrom data.table as.data.table
#' @export
ATL08_seg_attributes_dt_LAS <- function(atl08_dt, output) {
  longitude <- latitude <- h_canopy <- NA

  dt <- data.table::as.data.table(atl08_dt[, list(
    X = longitude,
    Y = latitude,
    Z = h_canopy
  )])

  dt_to_las(dt, output)
}
