#' Clips ICESat-2 ATL03 data
#'
#' @param atl03 [`icesat2.atl03_h5-class`] object, obtained through [`ATL03_read()`]
#' for clipping
#' @param bbox [`numeric-class`] or [`terra::SpatExtent`] for clipping.
#'
#' @return Returns the clipped S4 object of class [`icesat2.atl03_h5-class`]
#'
#' @description This function clips the ATl03 HDF5 file. This function
#' will only clip the beam groups within hdf5, it won't change metadata
#' or ancillary data.
#'
#' @examples
#' # Specifying the path to ICESat-2 ATL03 data (zip file)
#' outdir <- tempdir()
#' atl03_fp_zip <- system.file("extdata",
#'   "ATL0302_A_2019103030338_O01964_T05337_02_001_01_sub.zip",
#'   package = "rICESat2Veg"
#' )
#'
#' # Unzipping ICESat-2 ATL03 data
#' atl03_path <- unzip(atl03_fp_zip, exdir = outdir)
#'
#' # Reading ICESat-2 ATL03 data (h5 file)
#' atl03 <- ATL03_read(atl03_path = atl03_path)
#'
#' close(atl03)
#' @import hdf5r
#' @export
setGeneric(
  "ATL03_h5_clipBox",
  function(atl03, bbox) {
    standardGeneric("ATL03_h5_clipBox")
  }
)

#' @include class.icesat2.R
#' @importClassesFrom terra SpatExtent
setMethod(
  "ATL03_h5_clipBox",
  signature = c("icesat2.atl03_h5", "SpatExtent"),
  function(atl03, bbox) {
    print("clipping by ext")
  }
)

setMethod(
  "ATL03_h5_clipBox",
  signature = c("icesat2.atl03_h5", "numeric"),
  function(atl03, bbox) {
    print("clipping by bbox")
    bbox_ext <- terra::ext(bbox[c(2, 4, 3, 1)])
    ATL03_h5_clipBox(atl03, bbox_ext)
  }
)

getPhotonsMask <- function(atl03, bbox) {

}