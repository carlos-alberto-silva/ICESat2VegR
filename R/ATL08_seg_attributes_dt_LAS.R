#' Export ICESat-2 ATL08 Segment Data to LAS Files
#'
#' @description
#' Converts ATL08 segment-level vegetation attributes extracted with
#' \code{ATL08_seg_attributes_dt()} into one or more LAS files.
#'
#' @param atl08_dt An S4 object of class
#' \code{ICESat2VegR::icesat2.atl08_dt-class}, typically the output of
#' \code{ATL08_seg_attributes_dt()}.
#' @param output Character. Output LAS file path. The function creates one LAS
#' file per UTM zone in the WGS84 datum.
#'
#' @return Invisibly returns a character vector with the written LAS file paths.
#'
#' @details
#' This function exports ATL08 segment-level observations as LAS point clouds
#' using longitude, latitude, and canopy height \code{h_canopy}.
#'
#' Because LAS files require projected coordinates, the input geographic
#' coordinates are automatically reprojected to the appropriate UTM coordinate
#' reference system before export. If the data span multiple UTM zones, one LAS
#' file is generated for each zone.
#'
#' The internal UTM-zone assignment follows helper routines based on the
#' latitude/longitude to UTM conversion algorithms developed by Chuck Gantz.
#'
#'
#' @references
#' Gantz, C. Latitude/Longitude to UTM Conversion Algorithms.
#' \url{https://oceancolor.gsfc.nasa.gov/docs/ocssw/LatLong-UTMconversion_8cpp_source.html}
#'
#' @seealso
#' \code{\link{ATL08_seg_attributes_dt}},
#' \code{\link{dt_to_las}}
#'
#' @examples
#' \dontrun{
#' atl08_path <- system.file("extdata",
#'   "atl08_clip.h5",
#'   package = "ICESat2VegR"
#' )
#'
#' atl08_h5 <- ATL08_read(atl08_path = atl08_path)
#'
#' atl08_dt <- ATL08_seg_attributes_dt(atl08_h5)
#'
#' outputLas <- tempfile(fileext = ".las")
#'
#' ATL08_seg_attributes_dt_LAS(
#'   atl08_dt,
#'   outputLas
#' )
#'
#' close(atl08_h5)
#' }
#'
#' @include lasTools.R
#' @importFrom data.table as.data.table
#'
#' @export
ATL08_seg_attributes_dt_LAS <- function(atl08_dt, output) {

  longitude <- latitude <- h_canopy <- NULL

  atl08_dt <- data.table::as.data.table(atl08_dt)

  required_cols <- c("longitude", "latitude", "h_canopy")

  if (!all(required_cols %in% names(atl08_dt))) {
    stop(
      "atl08_dt must contain the columns: ",
      paste(required_cols, collapse = ", "),
      "."
    )
  }

  dt <- data.table::as.data.table(atl08_dt[, list(
    X = longitude,
    Y = latitude,
    Z = h_canopy
  )])

  dt_to_las(dt, output)
}
