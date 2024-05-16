#' Converts ATL03 photon cloud to LAS
#'
#' @param atl03_dt  An S4 object of class [`ICESat2VegR::icesat2.atl03_dt-class`]
#' (output of [`ATL03_photons_attributes_dt()`] function).
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
#'
#' # ATL03 file path
#' atl03_path <- system.file("extdata",
#'   "atl03_clip.h5",
#'   package = "ICESat2VegR"
#' )
#'
#' # Reading ATL03 data (h5 file)
#' atl03_h5 <- ATL03_read(atl03_path = atl03_path)
#'
#' # Extracting ATL03 and ATL08 photons and heights
#' atl03_dt <- ATL03_photons_attributes_dt(atl03_h5, beam = "gt1r")
#'
#' outdir <- tempdir()
#' ATL03_photons_attributes_dt_LAS(
#'   atl03_dt,
#'   file.path(outdir, "output.laz")
#' )
#'
#' close(atl03_h5)
#' @include utmTools.R lasTools.R
#' @importFrom data.table as.data.table
#' @export
ATL03_photons_attributes_dt_LAS <- function(atl03_dt, output) {
  # Define placeholders to avoid R CMD check warnings
  lon_ph <- lat_ph <- h_ph <- NA

  # Convert the input ATL03 data table to a data.table with required columns
  dt <- data.table::as.data.table(atl03_dt[, list(
    X = lon_ph,
    Y = lat_ph,
    Z = h_ph
  )])

  # Call the function to convert the data.table to LAS format and save to output path
  dt_to_las(dt, output)
}
