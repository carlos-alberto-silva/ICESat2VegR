#' Clip ICESat-2 ATL03 HDF5 Data Using Geometry-Based Boundaries
#'
#' @description
#' Clips an ICESat-2 ATL03 HDF5 file using one or more geometry-based clipping
#' objects provided as a [`terra::SpatVector`]. Each geometry in `vect`
#' defines an individual clipping region. The function iteratively clips the
#' ATL03 data for each region and writes a separate output HDF5 file for every
#' geometry. The name of each output file is automatically appended with the
#' value of the attribute specified in `split_by`.
#'
#' Only datasets within the ATL03 beam groups are spatially clipped; metadata
#' and ancillary non-beam groups remain unchanged in the resulting files.
#'
#' @param atl03 An [`ICESat2VegR::icesat2.atl03_h5-class`] object obtained using
#'   [`ATL03_read()`], representing the ATL03 HDF5 file to be clipped.
#'
#' @param output Character. Path to the output filename. The final written files
#'   will append the unique value of `split_by` for each clipping geometry,
#'   for example:
#'   `output_id1.h5`, `output_id2.h5`, etc.
#'
#' @param clip_obj A [`terra::SpatVector`] object containing the geometric clip
#'   boundaries. Each feature (row) in the SpatVector defines an independent
#'   clipping region.
#'
#' @param split_by Character. The SpatVector attribute whose values are used to
#'   identify and name each clipping geometry. Defaults to `id`.
#'
#' @param beam Character vector specifying which ATL03 beams to include. The
#'   default includes all six beams:
#'   `c("gt1l", "gt2l", "gt3l", "gt1r", "gt2r", "gt3r")`.
#'
#' @param additional_groups Character vector of non-beam HDF5 groups that should
#'   be copied unchanged into each clipped file. Defaults to:
#'   `c("orbit_info")`.
#'
#' @return
#' Returns a list of clipped S4 objects of class
#' [`ICESat2VegR::icesat2.atl03_h5-class`], one for each clipping geometry
#' contained in `vect`.
#'
#' @examples
#' # ATL03 file path
#' atl03_path <- system.file("extdata",
#'   "atl03_clip.h5",
#'   package = "ICESat2VegR"
#' )
#'
#' # Read ATL03 file
#' atl03_h5 <- ATL03_read(atl03_path)
#'
#' # Output base filename
#' output <- tempfile(fileext = ".h5")
#'
#' # Load clipping geometries
#' vect_path <- system.file("extdata", "clip_geom.shp",
#'   package = "ICESat2VegR"
#' )
#' vect <- terra::vect(vect_path)
#'
#' # Clip ATL03 data using polygon geometries
#' atl03_clipped_list <- ATL03_h5_clipGeometry(
#'   atl03_h5,
#'   output,
#'   vect,
#'   split_by = "id"
#' )
#'
#' close(atl03_h5)
#'
#' @import hdf5r
#' @include clipTools.R
#' @export
#' @export
ATL03_h5_clipGeometry <- function(
    atl03, output, clip_obj, split_by = NULL, beam = c("gt1r", "gt2r", "gt3r", "gt1l", "gt2l", "gt3l"),
    additional_groups = c("orbit_info")) {
  geom <- terra::aggregate(clip_obj)

  ATL03_h5_clip(atl03, output, geom, ATL03_segments_mask_geometry, beam, additional_groups)
}
