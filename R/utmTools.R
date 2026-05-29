# Function to convert Lat/Long to UTM Zone
# Credits to Chuck Gantz- chuck.gantz@globalstar.com
# https://oceancolor.gsfc.nasa.gov/docs/ocssw/LatLong-UTMconversion_8cpp_source.html
#' Convert latitude and longitude to UTM zone number
#'
#' @description
#' Computes the Universal Transverse Mercator (UTM) zone number for a given
#' set of geographic coordinates. Handles the special Norway zone correction
#' (zone 32 for parts of western Norway) and the Svalbard special zones
#' (31, 33, 35, and 37) as defined by the UTM standard.
#'
#' @param Lat Numeric vector of latitude values in decimal degrees.
#'   Valid range is -90 to 90.
#' @param Long Numeric vector of longitude values in decimal degrees.
#'   Valid range is -180 to 180.
#'
#' @return An integer vector of UTM zone numbers (1–60), the same length
#'   as the input vectors.
#'
#' @references
#' Chuck Gantz (chuck.gantz@@globalstar.com).
#' <https://oceancolor.gsfc.nasa.gov/docs/ocssw/LatLong-UTMconversion_8cpp_source.html>
#'
#' @examples
#' # Standard location (New York, USA) — zone 18
#' latLongToUtmZone(Lat = 40.71, Long = -74.01)
#'
#' # Norway special case (Bergen) — returns zone 32 instead of 31
#' latLongToUtmZone(Lat = 60.39, Long = 5.32)
#'
#' # Svalbard special case — returns zone 33
#' latLongToUtmZone(Lat = 78.22, Long = 15.65)
#'
#' # Multiple coordinates at once
#' latLongToUtmZone(
#'   Lat  = c(40.71, 60.39, 78.22, -33.87),
#'   Long = c(-74.01, 5.32, 15.65, 151.21)
#' )
#'
#' @export
latLongToUtmZone <- function(Lat, Long) {
  ZoneNumber <- as.integer((floor((Long + 180) / 6) %% 60) + 1)

  mask32 <- Lat >= 56.0 & Lat < 64.0 & Long >= 3.0 & Long < 12.0
  ZoneNumber[mask32] <- 32

  #########################
  # Special Svalbard zones
  #########################
  maskSvalbard <- Lat >= 72.0 & Lat < 84.0 & Long >= 0 & Long <= 42

  # Longitude intervals
  # [0,9) = 31
  # [9, 21) = 33
  # [21,33) = 35
  # [33,42) = 37
  svalbardIntervals <- base::findInterval(
    Long[maskSvalbard],
    c(0, 9, 21, 33, 42)
  )
  ZoneNumber[maskSvalbard] <- svalbardIntervals * 2 + 29

  ZoneNumber
}

#' Group coordinates by UTM/EPSG zone and return row index masks
#'
#' @description
#' For a set of geographic coordinates, computes the EPSG code of the
#' corresponding UTM projected coordinate reference system for each point
#' and returns a \code{data.table} grouping the row indices (\code{.I})
#' by EPSG code. This is useful for batch-reprojecting subsets of a
#' dataset, where each subset shares the same UTM zone.
#'
#' EPSG codes follow the pattern \code{326XX} for the northern hemisphere
#' and \code{327XX} for the southern hemisphere, where \code{XX} is the
#' two-digit UTM zone number (e.g., EPSG:32618 = WGS84 / UTM zone 18N).
#'
#' @param Lat Numeric vector of latitude values in decimal degrees.
#'   Valid range is -90 to 90.
#' @param Long Numeric vector of longitude values in decimal degrees.
#'   Valid range is -180 to 180.
#'
#' @return A \code{data.table} with one row per unique EPSG code found
#'   among the input coordinates. Columns are:
#'   \describe{
#'     \item{epsg}{Integer. The EPSG code of the UTM zone
#'       (e.g., 32618 for UTM zone 18N, 32718 for UTM zone 18S).}
#'     \item{mask}{List. Each element is an integer vector of row indices
#'       (positions in the original \code{Lat}/\code{Long} input) that
#'       fall within that EPSG zone.}
#'   }
#'
#' @examples
#' # Mix of northern and southern hemisphere points spanning two UTM zones
#' lats  <- c(40.71, 48.85, -33.87, -23.55)
#' longs <- c(-74.01, 2.35, 151.21, -46.63)
#'
#' result <- latLongToUtmMask(Lat = lats, Long = longs)
#' print(result)
#'
#' # Inspect which original rows belong to each zone
#' result$epsg    # EPSG codes
#' result$mask    # row indices per zone
#'
#' # Practical use: apply a different projection to each group
#' for (i in seq_len(nrow(result))) {
#'   idx  <- result$mask[[i]]
#'   zone_epsg <- result$epsg[i]
#'   message("EPSG ", zone_epsg, " covers rows: ", paste(idx, collapse = ", "))
#' }
#'
#' @import data.table
#' @export
latLongToUtmMask <- function(Lat, Long) {
  utmZones <- latLongToUtmZone(Lat, Long)
  isSouth <- Lat < 0
  epsg <- (326 + isSouth) * 100 + utmZones
  dt <- data.table::data.table(Lat, Long, epsg)
  dt <- dt[, list(mask = list(.I)), by = epsg]
  dt
}
