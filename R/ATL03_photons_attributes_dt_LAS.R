#' Converts ATL03 photon cloud to LAS
#'
#' @param atl03_dt  An S4 object of class [ICESat2VegR::icesat2.atl03_dt]
#' (output of [`ATL03_photons_attributes_dt()`] function).
#' @param output character. The output path of for the LAS(Z) file(s)
#' The function will create one LAS file per UTM Zone in WGS84 datum.
#'
#' @return Nothing, it just saves outputs as LAS file in disk
#'
#' @examples
#'
#' # Specifying the path to ATL03 and ATL08 file (zip file)
#' outdir <- tempdir()
#' atl03_zip <- system.file("extdata",
#'   "ATL03_20220401221822_01501506_005_01.zip",
#'   package = "ICESat2VegR"
#' )
#'
#' # Unzipping ATL03 file
#' atl03_path <- unzip(atl03_zip, exdir = outdir)
#'
#' # Reading ATL03 data (h5 file)
#' atl03_h5 <- ATL03_read(atl03_path = atl03_path)
#'
#' # # Extracting ATL03 and ATL08 photons and heights
#' atl03_dt <- ATL03_photons_attributes_dt(atl03_h5, beam = "gt1r")
#'
#' ATL03_photons_attributes_dt_LAS(
#'   atl03_dt,
#'   file.path(outdir, "output.laz")
#' )
#'
#' close(atl03_h5)
#' @include utmTools.R
#' @importFrom data.table as.data.table
#' @export
ATL03_photons_attributes_dt_LAS <- function(atl03_dt, output) {
  lon_ph <- lat_ph <- ph_h <- mask <- NA

  dt <- data.table::as.data.table(atl03_dt[, list(
    X = lon_ph,
    Y = lat_ph,
    Z =  h_ph
  )])
  # class(dt2)<-"data.table"

  names(dt) <- c("X", "Y", "Z")


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
