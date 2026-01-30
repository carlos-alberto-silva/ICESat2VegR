#' Read ICESat-2 ATL08 data
#'
#' @description
#' Read the ICESat-2 Land and Vegetation Along-Track Product (ATL08)
#' from an HDF5 (.h5) file.
#'
#' @param atl08_path
#'   File path to an ICESat-2 ATL08 dataset in HDF5 format (.h5).
#'
#' @return
#'   An S4 object of class [`icesat2.atl08_h5-class`] containing ICESat-2
#'   ATL08 data.
#'
#' @seealso
#'   \url{https://icesat-2.gsfc.nasa.gov/sites/default/files/page_files/ICESat2_ATL08_ATBD_r007.pdf}
#'
#' @examples
#' # Specify the path to an example ATL08 file
#' atl08_path <- system.file(
#'   "extdata",
#'   "atl08_clip.h5",
#'   package = "ICESat2VegR"
#' )
#'
#' # Read ICESat-2 ATL08 data (HDF5 file)
#' atl08 <- ATL08_read(atl08_path = atl08_path)
#' close(atl08)
#'
#' @import hdf5r
#' @include class.icesat2.R zzz.R
#' @export
setGeneric("ATL08_read", function(atl08_path) {
  standardGeneric("ATL08_read")
})


#' Read ICESat-2 ATL08 data from a local HDF5 file
#'
#' @description
#' Read the ICESat-2 Land and Vegetation Along-Track Product (ATL08)
#' from a local HDF5 (.h5) file.
#'
#' @param atl08_path
#'   Character. File path to ICESat-2 ATL08 data in HDF5 format (.h5).
#'
#' @return
#'   An S4 object of class [`icesat2.atl08_h5-class`] containing ICESat-2
#'   ATL08 data.
#'
#' @examples
#' # Specify the path to an example ATL08 file
#' atl08_path <- system.file(
#'   "extdata",
#'   "atl08_clip.h5",
#'   package = "ICESat2VegR"
#' )
#'
#' # Read ICESat-2 ATL08 data (HDF5 file)
#' atl08 <- ATL08_read(atl08_path = atl08_path)
#' close(atl08)
#'
#' @import hdf5r
#' @keywords internal
setMethod(
  "ATL08_read",
  signature = c("character"),
  function(atl08_path) {
    if (
      !(
        file.exists(atl08_path) &&
        tools::file_ext(atl08_path) == "h5"
      )
    ) {
      stop("atl08_path must be a path to a .h5 file")
    }

    atl08 <- ICESat2.h5_local$new(h5 = atl08_path)
    prepend_class(atl08, "icesat2.atl08_h5")
    return(atl08)
  }
)


#' Read ICESat-2 ATL08 data from multiple granules
#'
#' @description
#' This method stops execution when multiple granules are provided,
#' because the package currently supports only one granule at a time.
#'
#' @param atl08_path
#'   Object of class `icesat2.granules_cloud`. This parameter represents
#'   multiple ICESat-2 ATL08 granules.
#'
#' @return
#'   This method always stops with an error indicating that only one
#'   granule at a time is supported.
#'
#' @keywords internal
setMethod(
  "ATL08_read",
  signature = c("icesat2.granules_cloud"),
  function(atl08_path) {
    stop(
      "The package currently works with only one granule at a time.\n",
      "Please try again using a single granule, for example atl08_path[i]."
    )
  }
)


#' Read ICESat-2 ATL08 data from a single granule in the cloud
#'
#' @description
#' Read the ICESat-2 Land and Vegetation Along-Track Product (ATL08)
#' from a single granule accessible in the cloud.
#'
#' @param atl08_path
#'   Object of class `icesat2.granule_cloud`. This parameter represents
#'   a single ICESat-2 ATL08 granule in the cloud.
#'
#' @return
#'   An S4 object of class [`icesat2.atl08_h5-class`] containing ICESat-2
#'   ATL08 data.
#'
#' @keywords internal
setMethod(
  "ATL08_read",
  signature = c("icesat2.granule_cloud"),
  function(atl08_path) {
    atl08 <- ICESat2.h5_cloud$new(h5 = atl08_path)
    prepend_class(atl08, "icesat2.atl08_h5")
    return(atl08)
  }
)
