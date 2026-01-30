#' Clip ICESat-2 ATL08 HDF5 Data Using Geometry-Based Boundaries
#'
#' @description
#' Clips an ICESat-2 ATL08 HDF5 file using one or more geometry-based clipping
#' regions provided as a [`terra::SpatVector`]. Each feature (row) in
#' \code{clip_obj} defines a separate clipping boundary. For every geometry, a new
#' clipped ATL08 HDF5 file is created with the beam-level datasets spatially
#' filtered to include only segments within the clipping region.
#'
#' The resulting output filenames are automatically suffixed with the value of
#' the attribute specified in \code{split_by}, enabling batch generation of
#' multiple clipped ATL08 files from a single input.
#'
#' Metadata, ancillary groups, and orbit information are preserved unchanged in
#' every output file.
#'
#' @param ATL08 An [`ICESat2VegR::icesat2.ATL08_h5-class`] object obtained via
#'   [`ATL08_read()`], representing the ATL08 HDF5 file to be clipped.
#'
#' @param output Character string defining the base output filepath. The final
#'   output files will be named using:
#'   \code{paste0(output, "_", <split_by_value>, ".h5")}.
#'
#' @param clip_obj A [`terra::SpatVector`] containing the polygon or multipolygon
#'   geometries used as clipping regions. Each feature corresponds to one
#'   clipping output.
#'
#' @param split_by Character. Name of the attribute column in \code{clip_obj} used
#'   to uniquely identify each clipping region and to suffix output filenames.
#'   Defaults to \code{"id"}.
#'
#' @param beam Character vector specifying which ATL08 beams to include.
#'   Defaults to:
#'   \code{c("gt1l", "gt2l", "gt3l", "gt1r", "gt2r", "gt3r")}.
#'
#' @param additional_groups Character vector of additional non-beam groups to
#'   be copied unchanged into each output file. The default is:
#'   \code{c("orbit_info")}.
#'
#' @return
#' A list of clipped S4 objects of class
#' [`ICESat2VegR::icesat2.ATL08_h5-class`], one for each geometry in
#' \code{clip_obj}.
#'
#' @details
#' Only the beam-level groups (e.g., \code{gt1l}, \code{gt2r}, …) are spatially
#' filtered. All metadata, orbit information, and ancillary ATL08 groups are
#' preserved without modification.
#'
#' The function:
#' \enumerate{
#'   \item Computes polygon-based masks for each beam–segment pair.
#'   \item Removes all ATL08 segment datasets outside the clipping geometry.
#'   \item Writes each clipped subset into its own ATL08 HDF5 file.
#'   \item Returns the loaded results as S4 objects.
#' }
#'
#' This is the geometry-based companion to bounding-box clipping functions in
#' the package.
#'
#' @examples
#' # ATL08 file path
#' ATL08_path <- system.file("extdata",
#'   "ATL08_clip.h5",
#'   package = "ICESat2VegR"
#' )
#'
#' # Read ATL08 file
#' ATL08_h5 <- ATL08_read(ATL08_path)
#'
#' # Base output name (actual filenames will append split_by values)
#' output <- tempfile(fileext = ".h5")
#'
#' # Load clipping polygons
#' clip_obj_path <- system.file("extdata", "clip_objs.shp",
#'   package = "ICESat2VegR"
#' )
#' clip_obj <- terra::vect(vect_path)
#'
#' # Clip ATL08 using polygon geometries
#' ATL08_clipped_list <- ATL08_h5_clipGeometry(
#'   ATL08_h5,
#'   output,
#'   vect,
#'   split_by = "id"
#' )
#'
#' close(ATL08_h5)
#'
#' @import hdf5r
#' @include clipTools.R
#' @export
ATL08_h5_clipGeometry <- function(
    ATL08, output, clip_obj, split_by = NULL, beam = c("gt1r", "gt2r", "gt3r", "gt1l", "gt2l", "gt3l"),
    additional_groups = c("orbit_info")) {
  geom <- terra::aggregate(clip_obj)
plot(geom)
  ATL08_h5_clip(ATL08, output, geom, ATL08_segments_mask_geometry, beam, additional_groups)
}
