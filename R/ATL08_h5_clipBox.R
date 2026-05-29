#' Clip ICESat-2 ATL08 HDF5 data by bounding box
#'
#' @description
#' Clips an ATL08 HDF5 file to a given bounding extent, preserving
#' the original file structure, metadata, and ancillary data unchanged.
#' Returns the clipped data as a new \code{icesat2.atl08_h5} object.
#'
#' @param atl08 An object of class \code{icesat2.atl08_h5},
#'   typically produced by [ATL08_read()].
#' @param output Character. Full path to the output HDF5 file
#'   (must include the \code{.h5} extension).
#' @param clip_obj Bounding extent for clipping. Supported inputs:
#'   \itemize{
#'     \item A numeric vector of length 4:
#'           \code{c(ymax, xmin, ymin, xmax)} in decimal degrees.
#'     \item A \code{\link[terra]{SpatExtent}} object
#'           (\code{xmin, xmax, ymin, ymax}).
#'   }
#' @param beam Character vector of beam names to include.
#'   Default is all six beams:
#'   \code{c("gt1r", "gt2r", "gt3r", "gt1l", "gt2l", "gt3l")}.
#' @param additional_groups Character vector of non-beam HDF5 groups
#'   to copy into the output file. Default is \code{"orbit_info"}.
#'
#' @return An object of class \code{icesat2.atl08_h5} pointing to
#'   the clipped output HDF5 file.
#'
#' @details
#' The function clips the ATL08 HDF5 file beam by beam, retaining only
#' land segments whose coordinates fall within the specified bounding
#' extent. All beam datasets (photon-level and segment-level) are clipped
#' consistently. Non-beam groups (e.g. \code{orbit_info}) are copied
#' unchanged. The clipped file is written to \code{output} and returned
#' as a new \code{icesat2.atl08_h5} object.
#'
#' Note: \code{output} must be a full file path ending in \code{.h5}.
#' If the file already exists and is open, the function will fail with
#' an HDF5 error — always \code{close()} any open handles before
#' rerunning.
#'
#' @examples
#' # Specifying the path to ATL08 file
#' atl08_path <- system.file("extdata",
#'   "atl08_clip.h5",
#'   package = "ICESat2VegR"
#' )
#'
#' # Reading ATL08 data (h5 file)
#' atl08_h5 <- ATL08_read(atl08_path = atl08_path)
#'
#' # Output file path in temp directory
#' output <- file.path(tempdir(), "atl08_clipped.h5")
#'
#' # Clip using a terra SpatExtent
#' library(terra)
#' ext_obj <- terra::ext(-106.5709, -106.5699, 41.5315, 41.535)
#'
#' atl08_clipped <- ATL08_h5_clipBox(
#'   atl08_h5,
#'   output = output,
#'   clip_obj = ext_obj
#' )
#' atl08_clipped
#'
#' # Inspect the clipped segments
#' atl08_seg_clip <- ATL08_seg_attributes_dt(atl08_clipped)
#' head(atl08_seg_clip)
#' nrow(atl08_seg_clip)
#'
#' close(atl08_h5)
#' close(atl08_clipped)
#'
#' @import data.table hdf5r
#' @export
ATL08_h5_clipBox <- function(
    atl08, output, clip_obj,
    beam = c("gt1r", "gt2r", "gt3r", "gt1l", "gt2l", "gt3l"),
    additional_groups = c("orbit_info")) {
  if (inherits(clip_obj, "numeric")) {
    clip_obj <- terra::ext(clip_obj[2], clip_obj[4], clip_obj[3], clip_obj[1])
  }

  ATL08_h5_clip(atl08, output, clip_obj, landsegmentsMask_clip_obj, beam, additional_groups)
}

#' Clips ICESat-2 ATL08 data
#'
#' @param atl08 [`ICESat2VegR::icesat2.atl08_h5-class`] object, obtained through [`ATL08_read()`]
#' for clipping
#' @param output character. Path to the output h5 file.
#' @param clip_obj [`terra::SpatVector-class`] for clipping
#' @param split_by [`character-class`]. The attribute name used for identifying
#' the different clip_objs. Default is "id"
#' @param beam [`character-class`]. The vector of beams to include, default
#' all c("gt1l", "gt2l", "gt3l", "gt1r", "gt2r", "gt3r")
#' @param additional_groups [`character-class`]. Other addional groups that should be included, default
#' c("orbit_info")
#'
#' @return Returns a list of clipped S4 object of class [`ICESat2VegR::icesat2.atl08_h5-class`]
#'
#' @description This function clips ATL08 HDF5 file within beam groups,
#' but keeps metada and ancillary data the same.
#'
#' @examples
#' # ATL08 file path
#' atl08_path <- system.file("extdata",
#'   "atl08_clip.h5",
#'   package = "ICESat2VegR"
#' )
#'
#' # Reading ATL08 data (h5 file)
#' atl08_h5 <- ATL08_read(atl08_path = atl08_path)
#'
#' output <- tempfile(fileext = ".h5")
#'
#' vect_path <- system.file("extdata",
#'   "clip_geom.shp",
#'   package = "ICESat2VegR"
#' )
#'
#' vect <- terra::vect(vect_path)
#'
#' # Clipping ATL08 photons by boundary box extent
#' atl08_photons_dt_clip <- ATL08_h5_clipGeometry(
#'   atl08_h5,
#'   output,
#'   vect,
#'   split_by = "id"
#' )
#'
#' close(atl08_h5)
#' @import hdf5r
#' @export
ATL08_h5_clipGeometry <- function(
    atl08, output, clip_obj, split_by = "id",
    beam = c("gt1r", "gt2r", "gt3r", "gt1l", "gt2l", "gt3l"),
    additional_groups = c("orbit_info")) {
  geom <- terra::aggregate(clip_obj)
  ATL08_h5_clip(atl08, output, clip_obj = geom, landSegmentsMask_fn = landsegmentsMask_geom, beam, additional_groups)
}

landSegments_clip_obj <- function(beam, clip_obj) {
  latitude <- longitude <- NA
  .I <- data.table::.I

  land_segments_dt <- data.table::data.table(
    latitude = beam[["land_segments/latitude"]][],
    longitude = beam[["land_segments/longitude"]][]
  )

  landSegmentsDt <- land_segments_dt[, list(latitude, longitude, .I)][
    longitude >= clip_obj$xmin &
      longitude <= clip_obj$xmax &
      latitude >= clip_obj$ymin &
      latitude <= clip_obj$ymax
  ]

  landSegmentsDt
}

landsegmentsMask_clip_obj <- function(beam, clip_obj) {
  landSegmentsclip_obj <- landSegments_clip_obj(beam, clip_obj)

  landSegmentsclip_obj$I
}

landsegmentsMask_geom <- function(beam, geom) {
  clip_obj <- terra::ext(geom)
  landSegmentsclip_obj <- landSegments_clip_obj(beam, clip_obj)

  pts <- terra::vect(
    landSegmentsclip_obj,
    geom = c("longitude", "latitude"),
    crs = "epsg:4326"
  )

  res <- terra::intersect(pts, geom)
  res$I
}
