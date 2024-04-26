#' Clips ICESat-2 ATL03 H5 data
#'
#' @param atl03 [`ICESat2VegR::icesat2.atl03_h5-class`] object,
#' obtained through [`ATL03_read()`] for clipping
#' @param output character. Path to the output h5 file.
#' @param bbox [`numeric-class`] or [`terra::SpatExtent`] for clipping, the
#' order of the bbox is the default from NASA's ICESat-2 CMS searching:
#' (ul_lat, ul_lon, lr_lat, lr_lon).
#' @param beam [`character-class`]. The vector of beams to include, default
#' all c("gt1l", "gt2l", "gt3l", "gt1r", "gt2r", "gt3r")
#' @param additional_groups [`character-class`]. Other addional groups that should be included, default
#' c("METADATA", "orbit_info", "quality_assessment", "atlas_impulse_response", "ancillary_data")

#'
#' @return Returns the clipped S4 object of class [`ICESat2VegR::icesat2.atl03_h5-class`]
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
#' # Bounding rectangle coordinates
#' xmin <- -107.7
#' xmax <- -106.5
#' ymin <- 32.75
#' ymax <- 42.75
#'
#' # Clipping ATL03 photons  by boundary box extent
#' atl03_photons_dt_clip <- ATL03_h5_clipBox(atl03_h5, outdir, xmin, xmax, ymin, ymax)
#'
#' close(atl03_h5)
#' @import hdf5r
#' @include clipTools.R
#' @export
ATL03_h5_clipBox <- function(
    atl03, output, bbox, beam = c("gt1r", "gt2r", "gt3r", "gt1l", "gt2l", "gt3l"),
    additional_groups = c("METADATA", "orbit_info", "quality_assessment", "atlas_impulse_response", "ancillary_data")) {
  ATL03_h5_clip(atl03, output, bbox, ATL03_segments_mask, beam, additional_groups)
}
