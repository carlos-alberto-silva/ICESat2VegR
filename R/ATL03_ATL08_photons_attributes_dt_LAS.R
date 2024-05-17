#' Converts ATL03/ATL08 classified photon cloud to LAS
#'
#' @param atl03_atl08_dt  An S4 object of class [`ICESat2VegR::icesat2.atl08_dt-class`]
#'                        containing ATL03 and ATL08 data (output of
#'                        [ATL03_ATL08_photons_attributes_dt_join()] function).
#' @param output          character. The output path of for the LAS(Z) file(s). The
#'                        function will create one LAS file per UTM Zone in WGS84 datum.
#' @param normalized      logical, default TRUE. Whether the output should be
#'                        normalized LAS or raw altitude.
#'
#' @return Nothing, it just saves outputs as LAS file in disk
#'
#' @examples
#' outdir <- tempdir()
#'
#' # ATL03 file path
#' atl03_path <- system.file("extdata",
#'   "atl03_clip.h5",
#'   package = "ICESat2VegR"
#' )
#'
#' # ATL08 file path
#' atl08_path <- system.file("extdata",
#'   "atl08_clip.h5",
#'   package = "ICESat2VegR"
#' )
#'
#' # Reading ATL03 data (h5 file)
#' atl03_h5 <- ATL03_read(atl03_path = atl03_path)
#'
#' # Reading ATL08 data (h5 file)
#' atl08_h5 <- ATL08_read(atl08_path = atl08_path)
#'
#' # Extracting ATL03 and ATL08 photons and heights
#' atl03_atl08_dt <- ATL03_ATL08_photons_attributes_dt_join(atl03_h5, atl08_h5)
#'
#' ATL03_ATL08_photons_attributes_dt_LAS(
#'   atl03_atl08_dt,
#'   output = file.path(outdir, "output.laz"),
#'   normalized = TRUE
#' )
#'
#' close(atl03_h5)
#' close(atl08_h5)
#' @include utmTools.R
#' @importFrom data.table as.data.table
#' @export
ATL03_ATL08_photons_attributes_dt_LAS <- function(atl03_atl08_dt, output, normalized = TRUE) {
  stopifnot("lidR is not available" = require(lidR))

  lon_ph <- lat_ph <- classed_pc_flag <- mask <- NA

  atl03_atl08_dt <- na.omit(atl03_atl08_dt)

  dt <- data.table::as.data.table(atl03_atl08_dt[, list(
    X = lon_ph,
    Y = lat_ph,
    Z = get(c("h_ph", "ph_h")[normalized + 1]),
    Classification = as.integer(classed_pc_flag + 1)
  )])

  maskZones <- latLongToUtmMask(dt$Y, dt$X)

  message("====================================================")
  message(sprintf("The provided data will be splitted into %s UTM zones", nrow(maskZones)))
  message("====================================================")

  for (ii in seq_along(maskZones$epsg)) {
    dtLocal <- dt[maskZones[ii, mask][[1]]]
    epsgCode <- maskZones[ii, epsg]
    epsg <- sprintf("epsg:%s", epsgCode)
    coordinates <- terra::project(as.matrix(dtLocal[, list(X, Y)]), from = "epsg:4326", to = epsg)
    dtLocal$X <- coordinates[, 1]
    dtLocal$Y <- coordinates[, 2]

    header <- suppressWarnings(lidR::LASheader(dtLocal))
    lidR::epsg(header) <- epsgCode
    header@PHB[["X scale factor"]] <- 0.000001
    header@PHB[["Y scale factor"]] <- 0.000001
    header@PHB[["Z scale factor"]] <- 0.01
    las <- suppressWarnings(lidR::LAS(dtLocal, header))

    localOutput <- gsub("(\\.la[sz])", sprintf("_%s\\1", epsgCode), output)
    lidR::writeLAS(las, localOutput)
    message(sprintf("EPSG: %s - saved as %s", epsgCode, localOutput))
  }
}
