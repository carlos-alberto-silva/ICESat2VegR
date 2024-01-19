# Function to convert Lat/Long to UTM Zone
# Credits to Chuck Gantz- chuck.gantz@globalstar.com
# https://oceancolor.gsfc.nasa.gov/docs/ocssw/LatLong-UTMconversion_8cpp_source.html
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

latLongToUtmMask <- function(Lat, Long) {
  utmZones <- latLongToUtmZone(Lat, Long)
  isSouth <- Lat < 0
  epsg <- (326 + isSouth) * 100 + utmZones
  dt <- data.table::data.table(Lat, Long, epsg)
  dt <- dt[, list(mask = list(.I)), by = epsg]
  dt
}
