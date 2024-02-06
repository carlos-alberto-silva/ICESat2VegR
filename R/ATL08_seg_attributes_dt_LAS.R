#' Converts ATL08 segments to LAS
#'
#' @param atl08_dt  An S4 object of class [ICESat2VegR::icesat2.atl08_dt] containing ATL03 and ATL08 data
#' (output of [ATL03_ATL08_photons_attributes_dt_join()] function).
#' @param output character. The output path of for the LAS(Z) file(s)
#' The function will create one LAS file per UTM Zone in WGS84 datum.
#'
#' @return Nothing, it just saves outputs as LAS file in disk
#'
#' @examples
#'
#' # Specifying the path to ATL03 and ATL08 file (zip file)
#' outdir <- tempdir()
#'
#' atl08_zip <- system.file("extdata",
#'   "ATL08_20220401221822_01501506_005_01.zip",
#'   package = "ICESat2VegR"
#' )
#'
#' # Unzipping ATL08 file
#' atl08_path <- unzip(atl08_zip, exdir = outdir)
#'
#' # Reading ATL08 data (h5 file)
#' atl08_h5 <- ATL08_read(atl08_path = atl08_path)
#'
#' # # Extracting ATL03 and ATL08 photons and heights
#' atl08_dt <- ATL08_seg_attributes_dt(atl08_h5)
#'
#' ATL08_seg_attributes_dt_LAS(
#'   atl08_dt,
#'   file.path(outdir, "output.laz")
#' )
#'
#' close(atl08_h5)
#' @include utmTools.R
#' @importFrom data.table as.data.table
#' @export
ATL03_seg_attributes_dt_LAS <- function(atl08_dt, output) {
  longitude <- latitude <- h_canopy <- classed_pc_flag <- mask <- NA

  dt <- data.table::as.data.table(atl08_dt[, list(
    X = longitude,
    Y = latitude,
    Z = h_canopy
  )])

  maskZones <- latLongToUtmMask(dt$Y, dt$X)

  message("====================================================")
  message(sprintf("The provided data will be splitted into %s UTM zones", nrow(maskZones)))
  message("====================================================")

  for (ii in seq_along(maskZones)) {
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
    las <- suppressWarnings(lidR::LAS(dtLocal, header))

    localOutput <- gsub("(\\.la[sz])", sprintf("_%s\\1", epsgCode), output)
    message(sprintf("EPSG: %s - saved as %s", epsgCode, localOutput))
    lidR::writeLAS(las, localOutput)
  }
}
