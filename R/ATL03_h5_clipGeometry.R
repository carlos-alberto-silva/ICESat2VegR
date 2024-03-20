#' Clips ICESat-2 ATL03 H5 data
#'
#' @param atl03 [`icesat2.atl03_h5-class`] object,
#' obtained through [`ATL03_read()`] for clipping
#' @param output character. Path to the output h5 file, the attribute for polygons
#' will be appended to the file name.
#' @param vect [`terra::SpatVector-class`] for clipping
#' @param polygon_id [`character-class`]. The attribute name used for identifying
#' the different polygons. Default is "id"
#'
#' @return Returns a list of clipped S4 object of class [`icesat2.atl03_h5-class`]
#'
#' @description This function clips ATL03 HDF5 file within beam groups,
#' but keeps metada and ancillary data the same.
#'
#' @examples
##' # Specifying the path to ATL03 file (zip file)
#' outdir <- tempdir()
#' atl03_zip <- system.file("extdata",
#'   "atl03_20220401221822_01501506_005_01.zip",
#'   package = "ICESat2VegR"
#' )
#'
#' # Unzipping ATL03 file
#' atl03_path <- unzip(atl03_zip, exdir = outdir)
#'
#' # Reading ATL03 data (h5 file)
#' atl03_h5 <- atl03_read(atl03_path = atl03_path)
#'
#' output <- file.path(outdir, "clipped.h5")
#'
#' vect_path <- system.file("extdata",
#'   "polygons.shp",
#'   package = "ICESat2VegR"
#' )
#'
#' vect <- terra::vect()
#'
#' # Clipping ATL03 photons  by boundary box extent
#' atl03_photons_dt_clip <- ATL03_h5_clipGeometry(
#'   atl03_h5,
#'   output,
#'   vect,
#'   polygon_id = "id"
#' )
#'
#' close(atl03_h5)
#' @import hdf5r
#' @include clipTools.R
#' @export
ATL03_h5_clipGeometry <- function(atl03, output, vect, polygon_id = NULL) {
  geom <- terra::union(vect)

  ATL03_h5_clip(atl03, output, geom, ATL03_segments_mask_geometry)
}
