#' Clip ICESat-2 ATL03 HDF5 Data Using a Bounding Extent
#'
#' @description
#' Clips an ICESat-2 ATL03 HDF5 file to a specified spatial extent. Only
#' geolocated photon and segment datasets within the ATL03 beam groups are
#' clipped; all metadata, orbit information, and ancillary groups are preserved
#' unchanged. This function provides bounding-box-based clipping, complementing
#' geometry-based clipping workflows.
#'
#' @param atl03 An [`ICESat2VegR::icesat2.atl03_h5-class`] object created by
#'   [`ATL03_read()`], representing the ATL03 HDF5 file to be clipped.
#'
#' @param output Character. Path to the output HDF5 file that will store the
#'   clipped ATL03 dataset.
#'
#' @param clip_obj Bounding extent used for clipping. Supported inputs:
#'   \itemize{
#'     \item A numeric vector of length 4:
#'       \code{c(ul_lat, ul_lon, lr_lat, lr_lon)}, following the ICESat-2 CMS
#'       bounding-box convention (upper-left latitude/longitude and
#'       lower-right latitude/longitude).
#'     \item A [`terra::SpatExtent`] object.
#'   }
#'
#' @param beam Character vector specifying which ATL03 beams to include.
#'   Defaults to:
#'   \code{c("gt1l", "gt2l", "gt3l", "gt1r", "gt2r", "gt3r")}.
#'
#' @param additional_groups Character vector of additional non-beam HDF5 groups
#'   to copy unchanged into the output file. Defaults to:
#'   \code{c("orbit_info")}.
#'
#' @return
#' An S4 object of class [`ICESat2VegR::icesat2.atl03_h5-class`] representing
#' the clipped ATL03 HDF5 file saved at the \code{output} location.
#'
#' @details
#' The function performs spatial subsetting at the photon and segment level.
#' Internally, it:
#'
#' \enumerate{
#'   \item Copies file- and group-level attributes.
#'   \item Extracts beam-level datasets and evaluates segment-level and
#'         photon-level masks using the supplied bounding box.
#'   \item Applies these masks to all datasets whose dimensions correspond to
#'         segments or photons.
#'   \item Recomputes dependent indexing datasets (e.g., \code{ph_index_beg})
#'         to maintain ATL03 structural integrity.
#'   \item Writes the clipped data into a new, valid ATL03 HDF5 file.
#' }
#'
#' This function preserves the full ATL03 metadata, ancillary information, and
#' file structure while trimming the spatial domain to the user-specified
#' bounding extent.
#'
#' @examples
#' # ATL03 file path
#' atl03_path <- system.file("extdata",
#'   "atl03_clip.h5",
#'   package = "ICESat2VegR"
#' )
#'
#' # Read ATL03 data (HDF5)
#' atl03_h5 <- ATL03_read(atl03_path)
#'
#' # Bounding rectangle coordinates (ul_lat, ul_lon, lr_lat, lr_lon)
#' xmin <- -106.5723
#' xmax <- -106.5693
#' ymin <- 41.533
#' ymax <- 41.537
#'
#' bbox <- c(ymax, xmin, ymin, xmax)
#'
#' # Clip ATL03 using the bounding extent
#' output <- tempfile(fileext = ".h5")
#' atl03_clip <- ATL03_h5_clipBox(
#'   atl03_h5,
#'   output = output,
#'   clip_obj = bbox
#' )
#'
#' close(atl03_h5)
#'
#' @import hdf5r
#' @include clipTools.R
#' @export
ATL03_h5_clipBox <- function(
    atl03, output, clip_obj, beam = c("gt1r", "gt2r", "gt3r", "gt1l", "gt2l", "gt3l"),
    additional_groups = c("orbit_info")) {
  if (inherits(clip_obj, "numeric")) {
    clip_obj <- terra::ext(clip_obj[2], clip_obj[4], clip_obj[3], clip_obj[1])
  }

  ATL03_h5_clip(atl03, output, clip_obj, ATL03_segments_mask, beam, additional_groups)
}
