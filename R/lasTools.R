#' @include utmTools.R
dt_to_las <- function(dt, output) {
  stopifnot("lidR is not available" = require(lidR))
  mask <- NA

  names(dt) <- c("X", "Y", "Z")
  maskZones <- latLongToUtmMask(dt$Y, dt$X)

  message("====================================================")
  message(sprintf("The provided data will be splitted into %s UTM zones", nrow(maskZones)))
  message("====================================================")

  for (ii in seq_len(nrow(maskZones))) {
    dtLocal <- dt[maskZones[ii, mask][[1]]]
    epsgCode <- maskZones[ii, epsg]
    epsg <- sprintf("epsg:%s", epsgCode)
    coordinates <- terra::project(as.matrix(dtLocal[, list(X, Y)]), from = "epsg:4326", to = epsg)
    dtLocal$X <- coordinates[, 1]
    dtLocal$Y <- coordinates[, 2]

    header <- suppressWarnings(lidR::LASheader(dtLocal))
    lidR::epsg(header) <- epsgCode
    header@PHB[["X scale factor"]] <- 0.01
    header@PHB[["Y scale factor"]] <- 0.01
    header@PHB[["Z scale factor"]] <- 0.01
    las <- suppressWarnings(lidR::LAS(
      dtLocal,
      header,
      crs = lidR::st_crs(epsg)
    ))

    localOutput <- gsub("(\\.la[sz])", sprintf("_%s\\1", epsgCode), output)
    message(sprintf("EPSG: %s - saved as %s", epsgCode, localOutput))
    lidR::writeLAS(las, localOutput)
  }
}