#' Export ICESat-2 ATL08 photon attributes to LAS files
#'
#' @description
#' Converts ICESat-2 ATL08 photon-level data stored as a data.table to one or
#' more LAS files.
#'
#' @param atl08_dt A data.table containing ICESat-2 ATL08 photon attributes.
#' The table must contain longitude, latitude, and photon height columns.
#' Expected column names are \code{longitude}, \code{latitude}, and \code{ph_h}.
#' @param output Character. Output LAS file path. If the data span multiple UTM
#' zones, one LAS file is created per UTM zone.
#'
#' @return Invisibly returns a character vector with the written LAS file paths.
#'
#' @details
#' This function is designed for ATL08 photon-level data, including classified
#' photon attributes such as photon height, classification flag, and photon
#' segment ID. The input must contain longitude, latitude, and photon height
#' information. Because LAS files require projected coordinates, the input
#' longitude and latitude values are automatically split by UTM zone and
#' reprojected before writing the LAS file. The internal UTM-zone assignment
#' follows helper routines based on the latitude/longitude to UTM conversion
#' algorithms developed by Chuck Gantz.
#'
#' @references
#' Gantz, C. Latitude/Longitude to UTM Conversion Algorithms.
#' \url{https://oceancolor.gsfc.nasa.gov/docs/ocssw/LatLong-UTMconversion_8cpp_source.html}

#' @examples
#' \dontrun{
#' atl08_h5 <- ATL08_read(atl08_path)
#'
#' atl08_photons <- ATL08_photons_attributes_dt(atl08_h5)
#'
#' ATL08_photons_attributes_dt_LAS(
#'   atl08_dt = atl08_photons,
#'   output = "atl08_photons.las"
#' )
#'
#' close(atl08_h5)
#' }
#'
#' @include lasTools.R
#' @importFrom data.table as.data.table
#' @export
ATL08_photons_attributes_dt_LAS <- function(atl08_dt, output) {

  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required.")
  }

  atl08_dt <- data.table::as.data.table(atl08_dt)

  required_cols <- c("longitude", "latitude", "ph_h")

  if (!all(required_cols %in% names(atl08_dt))) {
    stop(
      "atl08_dt must contain the columns: longitude, latitude, and ph_h. ",
      "If your ATL08 photon table does not contain longitude and latitude, ",
      "join the ATL08 photons with ATL03 geolocation information first."
    )
  }

  longitude <- latitude <- ph_h <- NULL

  dt <- data.table::as.data.table(atl08_dt[, list(
    X = longitude,
    Y = latitude,
    Z = ph_h
  )])

  dt_to_las(dt, output)
}
