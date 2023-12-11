#' @importFrom hdf5r H5File
setRefClass("H5File")
requireNamespace("data.table")

#' Class for ICESat-2 ATL08
#'
#' @slot h5 Object of class [`H5File`][hdf5r::H5File-class] from `hdf5r` package containing the
#' ICESat-2 Land and Vegetation Along-Track Products (ATL08)

#'
#' @seealso [`H5File`][hdf5r::H5File-class] in the `hdf5r` package and
#' \url{https://icesat-2.gsfc.nasa.gov/sites/default/files/page_files/ICESat2_ATL08_ATBD_r006.pdf}
#'
#' @import methods
#' @export
icesat2.atl08_h5 <- setClass(
  Class = "icesat2.atl08_h5",
  slots = list(h5 = "H5File")
)

#' @importFrom hdf5r H5File
setRefClass("H5File")
requireNamespace("data.table")

#' Class for ICESat-2 ATL03
#'
#' @slot h5 Object of class [`H5File`][hdf5r::H5File-class] from `hdf5r` package containing the
#' ICESat-2 Global Geolocated Photon Data (ATL03)

#'
#' @seealso [`H5File`][hdf5r::H5File-class] in the `hdf5r` package and
#' \url{https://icesat-2.gsfc.nasa.gov/sites/default/files/page_files/ICESat2_ATL03_ATBD_r006.pdf}
#'
#' @import methods
#' @export
icesat2.atl03_h5 <- setClass(
  Class = "icesat2.atl03_h5",
  slots = list(h5 = "H5File")
)


#' @importFrom data.table data.table
setRefClass("data.table")
requireNamespace("data.table")

#' Class for ATL08 attributes
#'
#' @slot data.table Object of class [`data.table`][data.table::data.table-class] from `data.table` package containing the
#' Extracted ICESat-2 ATL08 attributes

#'
#' @seealso [`data.table`][data.table::data.table-class] in the `data.table` package and
#' \url{https://icesat-2.gsfc.nasa.gov/sites/default/files/page_files/ICESat2_ATL03_ATBD_r006.pdf}
#'
#' @import methods
#' @export
icesat2.atl08_dt <- setClass(
  Class = "icesat2.atl08_dt",
  slots = list(dt = "data.table")
)

#' @importFrom data.table data.table
setRefClass("data.table")
requireNamespace("data.table")

#' Class for joined ATL03 and ATL08 attributes
#'
#' @slot data.table Object of class [`data.table`][data.table::data.table-class] from `data.table` package containing the
#' Joined ICESat-2 ATL03 ATL08 attributes

#'
#' @seealso [`data.table`][data.table::data.table-class] in the `data.table` package and
#' \url{https://icesat-2.gsfc.nasa.gov/sites/default/files/page_files/ICESat2_ATL03_ATBD_r006.pdf}
#'
#' @import methods
#' @export
icesat2.atl03atl08_dt <- setClass(
  Class = "icesat2.atl03atl08_dt",
  slots = list(dt = "data.table")
)
