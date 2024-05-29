#' Read ICESat-2 ATL08 data
#'
#' @description This function reads the ICESat-2 Land and
#' Vegetation Along-Track Products (ATL08) as h5 file.
#'
#' @param atl08_path File path pointing to ICESat-2 ATL08 data. Data in HDF5 Hierarchical Data Format (.h5).
#'
#' @return Returns an S4 object of class [`icesat2.atl08_h5-class`] containing ICESat-2 ATL08 data.
#'
#' @seealso \url{https://icesat-2.gsfc.nasa.gov/sites/default/files/page_files/ICESat2_ATL08_ATBD_r006.pdf}
#'
#' @examples
#' # Specifying the path to ATL08 file
#' atl08_path <- system.file("extdata",
#'   "atl08_clip.h5",
#'   package = "ICESat2VegR"
#' )
#'
#' # Reading ICESat-2 ATL08 data (h5 file)
#' atl08 <- ATL08_read(atl08_path = atl08_path)
#' close(atl08)
#' @import hdf5r
#' @include class.icesat2.R zzz.R
#' @export
setGeneric("ATL08_read", function(atl08_path) {
  standardGeneric("ATL08_read")
})


#' Read ICESat-2 ATL08 data from local HDF5 file
#'
#' @description This method reads the ICESat-2 Land and
#' Vegetation Along-Track Products (ATL08) from a local HDF5 file.
#'
#' @param atl08_path Character. File path pointing to ICESat-2 ATL08 data in HDF5 format (.h5).
#'
#' @return Returns an S4 object of class [`icesat2.atl08_h5-class`] containing ICESat-2 ATL08 data.
#'
#' @examples
#' # Specifying the path to ATL08 file
#' atl08_path <- system.file("extdata",
#'   "atl08_clip.h5",
#'   package = "ICESat2VegR"
#' )
#'
#' # Reading ICESat-2 ATL08 data (h5 file)
#' atl08 <- ATL08_read(atl08_path = atl08_path)
#' close(atl08)
#' @import hdf5r
#' @export
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
      stop("atl08_path must be a path to a h5 file")
    }

    atl08 <- ICESat2.h5_local$new(h5 = atl08_path)
    prepend_class(atl08, "icesat2.atl08_h5")
    return(atl08)
  }
)

#' Read ICESat-2 ATL08 data from multiple granules
#'
#' @description This method stops the process when multiple granules are provided.
#'
#' @param atl08_path Object of class `icesat2.granules_cloud`. This parameter represents multiple granules of ICESat-2 ATL08 data.
#'
#' @return This method stops execution with an error message indicating that the package works with only one granule at a time.
#'
#' @examples
#' # Specifying the path to ATL08 file
#' atl08_path <- system.file("extdata",
#'   "atl08_clip.h5",
#'   package = "ICESat2VegR"
#' )
#'
#' # Reading ICESat-2 ATL08 data (h5 file)
#' atl08 <- ATL08_read(atl08_path = atl08_path)
#' close(atl08)
#' @export
setMethod(
  "ATL08_read",
  signature = c("icesat2.granules_cloud"),
  function(atl08_path) {
    stop("For now the package only works with one granule at a time
try with only one granule [i].")
  }
)

#' Read ICESat-2 ATL08 data from a single granule in the cloud
#'
#' @description This method reads the ICESat-2 Land and Vegetation Along-Track Products (ATL08) from a single granule in the cloud.
#'
#' @param atl08_path Object of class `icesat2.granule_cloud`. This parameter represents a single granule of ICESat-2 ATL08 data.
#'
#' @return Returns an S4 object of class [`icesat2.atl08_h5-class`] containing ICESat-2 ATL08 data.
#'
#' @examples
#' # Specifying the path to ATL08 file
#' atl08_path <- system.file("extdata",
#'   "atl08_clip.h5",
#'   package = "ICESat2VegR"
#' )
#'
#' # Reading ICESat-2 ATL08 data (h5 file)
#' atl08 <- ATL08_read(atl08_path = atl08_path)
#' close(atl08)
#' @export
setMethod(
  "ATL08_read",
  signature = c("icesat2.granule_cloud"),
  function(atl08_path) {
    atl08 <- ICESat2.h5_cloud$new(h5 = atl08_path)
    prepend_class(atl08, "icesat2.atl08_h5")
    return(atl08)
  }
)
