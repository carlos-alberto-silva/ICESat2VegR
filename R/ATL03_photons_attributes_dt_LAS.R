#' Export ICESat-2 ATL03 Photon Data to LAS Files
#'
#' @description
#' Converts ATL03 photon-level data extracted with
#' \code{ATL03_photons_attributes_dt()} into one or more LAS files.
#'
#' @param atl03_dt An S4 object of class
#' \code{ICESat2VegR::icesat2.atl03_dt-class}, usually the output of
#' \code{ATL03_photons_attributes_dt()}.
#' @param output Character. Output LAS file path. The function creates one LAS
#' file per UTM zone in the WGS84 datum.
#'
#' @return Invisibly returns a character vector with the written LAS file paths.
#'
#' @details
#' LAS files require projected coordinates. This function uses internal helper
#' functions to identify the appropriate UTM zone for each photon, reproject
#' the original ICESat-2 geographic coordinates, and write one LAS file per UTM
#' zone using the package internal LAS writer. The internal UTM-zone assignment
#' follows helper routines based on the latitude/longitude to UTM conversion
#' algorithms developed by Chuck Gantz.
#'
#' @references
#' Gantz, C. Latitude/Longitude to UTM Conversion Algorithms.
#' \url{https://oceancolor.gsfc.nasa.gov/docs/ocssw/LatLong-UTMconversion_8cpp_source.html}
#'
#' @examples
#' \dontrun{
#' atl03_path <- system.file("extdata",
#'   "atl03_clip.h5",
#'   package = "ICESat2VegR"
#' )
#'
#' atl03_h5 <- ATL03_read(atl03_path = atl03_path)
#'
#' atl03_dt <- ATL03_photons_attributes_dt(atl03_h5, beam = "gt1r")
#'
#' outdir <- tempdir()
#'
#' ATL03_photons_attributes_dt_LAS(
#'   atl03_dt,
#'   file.path(outdir, "atl03_photons.las")
#' )
#'
#' close(atl03_h5)
#' }
#'
#' @include utmTools.R lasTools.R
#' @importFrom data.table as.data.table
#'
#' @export
ATL03_photons_attributes_dt_LAS <- function(atl03_dt, output) {

  lon_ph <- lat_ph <- h_ph <- NULL

  dt <- data.table::as.data.table(atl03_dt[, list(
    X = lon_ph,
    Y = lat_ph,
    Z = h_ph
  )])

  dt_to_las(dt, output)
}
