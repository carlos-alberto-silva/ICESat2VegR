#' Read ICESat-2 ATL08 data
#'
#' @description This function reads the ICESat-2 Land and
#' Vegetation Along-Track Products (ATL08) as h5 file.
#'
#' @param atl08_path File path pointing to ICESat-2 ATL08 data. Data in HDF5 Hierarchical Data Format (.h5).
#'
#' @return Returns an S4 object of class [`icesat2.atl08_dt-class`] containing ICESat-2 ATL08 data.
#'
#' @seealso \url{https://icesat-2.gsfc.nasa.gov/sites/default/files/page_files/ICESat2_ATL08_ATBD_r006.pdf}
#'
#' @examples
#' # ATL08 file path
#' atl08_path <- system.file("extdata",
#'   "atl08_clip.h5",
#'   package = "ICESat2VegR"
#' )
#'
#' # Reading ICESat-2 ATL08 data (h5 file)
#' atl08 <- ATL08_read(atl08_path = atl08_path)
#'
#' close(atl08)
#' @import hdf5r
#' @export
ATL08_attributes_dt_LAS <- function(atl08_path) {
  if (!is.character(atl08_path) | !tools::file_ext(atl08_path) == "h5") {
    stop("atl08_path must be a path to a h5 file")
  }

  atl08_h5 <- hdf5r::H5File$new(atl08_path, mode = "r")
  atl08 <- new("icesat2.atl08", h5 = atl08_h5)
  return(atl08)
}
