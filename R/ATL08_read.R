#' Read ICESat-2 ATL08 data
#'
#' @description This function reads the ICESat-2 Land and
#' Vegetation Along-Track Products (ATL08) as h5 file.
#'
#'
#' @usage ATL08_read(atl08_path)
#'
#' @param atl08_path File path pointing to ICESat-2 ATL08 data. Data in HDF5 Hierarchical Data Format (.h5).
#'
#' @return Returns an S4 object of class ["icesat2.atl08_h5"] containing ICESat-2 ATL08 data.
#'
#' @seealso \url{https://icesat-2.gsfc.nasa.gov/sites/default/files/page_files/ICESat2_ATL08_ATBD_r006.pdf}
#'
#' @examples
#' # Specifying the path to ICESat-2 ATL08 data (zip file)
#' outdir <- tempdir()
#' atl08_fp_zip <- system.file("extdata",
#'   "ATL0802_A_2019108080338_O01964_T05337_02_001_01_sub.zip",
#'   package = "ICESat2VegR"
#' )
#'
#' # Unzipping ICESat-2 ATL08 data
#' atl08_path <- unzip(atl08_fp_zip, exdir = outdir)
#'
#' # Reading ICESat-2 ATL08 data (h5 file)
#' atl08 <- ATL08_read(atl08_path = atl08_path)
#'
#' close(atl08)
#' @import hdf5r
#' @export
setGeneric("ATL08_read", function(atl08_path) {
  standardGeneric("ATL08_read")
})


#' @include class.icesat2.R
#' @include zzz.R


setMethod(
  "ATL08_read",
  signature = c("ANY"),
  function(atl08_path) {
    if (
      !(
        file.exists(atl08_path) &&
          tools::file_ext(atl08_path) == "h5"
      )
    ) {
      stop("atl08_path must be a path to a h5 file")
    }

    atl08_h5 <- hdf5r::H5File$new(atl08_path, mode = "r")
    prepend_class(atl08_h5, "icesat2.hdf5r")

    atl08 <- new("icesat2.atl08_h5", h5 = atl08_h5)
    return(atl08)
  }
)

setMethod(
  "ATL08_read",
  signature = c("icesat2.granule_cloud"),
  function(atl08_path) {
    if (inherits(h5py, "logical")) {
      if (!reticulate::py_module_available("h5py")) {
        reticulate::py_install("h5py")
      }
      h5py <<- reticulate::import("h5py")
    }
    granules <- earthaccess$open(list(atl08_path))
    atl08_h5 <- h5py$File(reticulate::py_to_r(granules)[[1]])
    prepend_class(atl08_h5, "icesat2.h5_cloud")

    atl08 <- new("icesat2.atl08_h5", h5 = atl08_h5)

    return(atl08)
  }
)
