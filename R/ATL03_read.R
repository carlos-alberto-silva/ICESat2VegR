#' Read ICESat-2 ATL03 data
#'
#' @description This function reads the ICESat-2 Global Geolocated Photons (ATL03) Product (ATL03) as h5 file.
#'
#' @usage ATL03_read(atl03_path)
#'
#' @param atl03_path Either file path pointing to ICESat-2 ATL03 h5 data
#' or a granule resulting from [`ATLAS_dataFinder_cloud()`].
#'
#' @return Returns an S4 object of class [`icesat2.atl03_dt`] containing ICESat-2 ATL03 data.
#'
#' @seealso \url{https://icesat-2.gsfc.nasa.gov/sites/default/files/page_files/ICESat2_ATL03_ATBD_r006.pdf}
#'
#' @examples
#' # Specifying the path to ICESat-2 ATL03 data (zip file)
#' outdir <- tempdir()
#' ATL03_fp_zip <- system.file("extdata",
#'   "ATL0302_A_2019108080338_O01964_T05337_02_001_01_sub.zip",
#'   package = "ICESat2VegR"
#' )
#'
#' # Unzipping ICESat-2 ATL03 data
#' atl03_path <- unzip(ATL03_fp_zip, exdir = outdir)
#'
#' # Reading ICESat-2 ATL03 data (h5 file)
#' ATL03 <- ATL03_read(atl03_path = atl03_path)
#'
#' close(ATL03)
#' @import hdf5r
#' @export
setGeneric("ATL03_character", function(atl03_path) {
  standardGeneric("ATL03_read")
})


#' @include zzz.R
#' @include class.icesat2.R
setMethod(
  "ATL03_read",
  signature = c("ANY"),
  function(atl03_path) {
    if (
      !(
        file.exists(atl03_path) &&
          tools::file_ext(atl03_path) == "h5"
      )
    ) {
      stop("atl03_path must be a path to a h5 file")
    }

    atl03 <- ICESat2.h5_local$new(h5 = atl03_path)
    prepend_class(atl03, "icesat2.atl03_h5")
    return(atl03)
  }
)

setMethod(
  "ATL03_read",
  signature = c("icesat2.granule_cloud"),
  function(atl03_path) {
    atl03 <- ICESat2.h5_cloud$new(h5 = atl03_path)
    prepend_class(atl03, "icesat2.atl03_h5")
    return(atl03)
  }
)
