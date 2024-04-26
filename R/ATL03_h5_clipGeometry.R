#' Clips ICESat-2 ATL03 H5 data
#'
#' @param atl03 [`ICESat2VegR::icesat2.atl03_h5-class`] object,
#' obtained through [`ATL03_read()`] for clipping
#' @param output character. Path to the output h5 file, the attribute for polygons
#' will be appended to the file name.
#' @param vect [`terra::SpatVector-class`] for clipping
#' @param polygon_id [`character-class`]. The attribute name used for identifying
#' the different polygons. Default is "id"
#' @param beam [`character-class`]. The vector of beams to include, default
#' all c("gt1l", "gt2l", "gt3l", "gt1r", "gt2r", "gt3r")
#' @param additional_groups [`character-class`]. Other addional groups that should be included, default
#' c("METADATA", "orbit_info", "quality_assessment", "atlas_impulse_response", "ancillary_data")
#'
#' @return Returns a list of clipped S4 object of class [`ICESat2VegR::icesat2.atl03_h5-class`]
#'
#' @description This function clips ATL03 HDF5 file within beam groups,
#' but keeps metada and ancillary data the same.
#'
#' @examples
#' # ATL03 file path
#' atl03_path <- system.file("extdata",
#'   "atl03_clip.h5",
#'   package = "ICESat2VegR"
#' )
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
ATL03_h5_clipGeometry <- function(
    atl03, output, vect, polygon_id = NULL, beam = c("gt1r", "gt2r", "gt3r", "gt1l", "gt2l", "gt3l"),
    additional_groups = c("METADATA", "orbit_info", "quality_assessment", "atlas_impulse_response", "ancillary_data")) {
  geom <- terra::union(vect)

  ATL03_h5_clip(atl03, output, geom, ATL03_segments_mask_geometry, beam, additional_groups)
}
