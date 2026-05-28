#' Clip ICESat-2 ATL08 HDF5 Data Using a Bounding Extent
#'
#' @description
#' Clips an ICESat-2 ATL08 HDF5 file to a specified spatial extent. Only
#' segment-level datasets within the ATL08 beam groups are clipped; all metadata,
#' orbit information, and ancillary groups are preserved unchanged. This
#' function is the bounding-box-based counterpart to geometry-based clipping
#' functions.
#'
#' @param atl08 An [`ICESat2VegR::icesat2.atl08_h5-class`] object obtained via
#'   [`ATL08_read()`], representing the ATL08 HDF5 file to be clipped.
#'
#' @param output Character path specifying where the clipped HDF5 file should be
#'   written.
#'
#' @param clip_obj Bounding extent used for clipping. Supported inputs:
#'   \itemize{
#'     \item A numeric vector of length 4:
#'       `c(ul_lat, ul_lon, lr_lat, lr_lon)`, following the ICESat-2 CMS
#'       bounding-box convention (upper-left latitude/longitude and
#'       lower-right latitude/longitude).
#'     \item A [`terra::SpatExtent`] object.
#'   }
#'
#' @param beam Character vector specifying which ATL08 beams to include.
#'   Defaults to:
#'   `c("gt1l", "gt2l", "gt3l", "gt1r", "gt2r", "gt3r")`.
#'
#' @param additional_groups Character vector specifying additional non-beam
#'   HDF5 groups to copy unchanged into the output file. Defaults to:
#'   `c("orbit_info")`.
#'
#' @return
#' An S4 object of class [`ICESat2VegR::icesat2.atl08_h5-class`] representing
#' the clipped ATL08 HDF5 file.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Copies file-level attributes and all requested non-beam groups.
#'   \item Clips beam-level datasets based on latitude/longitude coordinates and
#'         the supplied spatial extent.
#'   \item Reconstructs dependent indexing datasets (e.g., photon index ranges)
#'         when necessary to maintain valid ATL08 structure.
#' }
#'
#' The resulting HDF5 file remains fully compliant with the ATL08 product
#' structure, but contains only data within the specified bounding extent.
#'
#' @examples
#' # ATL08 file path
#' atl08_path <- system.file("extdata",
#'   "atl08_clip.h5",
#'   package = "ICESat2VegR"
#' )
#'
#' # Read ATL08 data
#' atl08_h5 <- ATL08_read(atl08_path)
#'
#' # Bounding rectangle coordinates (ul_lat, ul_lon, lr_lat, lr_lon)
#' ul_lon <- -106.5723
#' lr_lon <- -106.5693
#' lr_lat <- 41.533
#' ul_lat <- 41.537
#'
#' # Clip ATL08 data using the bounding extent
#' atl08_clip <- ATL08_h5_clipBox(
#'   atl08_h5,
#'   output = tempfile(fileext = ".h5"),
#'   clip_obj = c(ul_lat, ul_lon, lr_lat, lr_lon)
#' )
#'
#' close(atl08_h5)
#' close(atl08_clip)
#'
#' @import hdf5r
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
