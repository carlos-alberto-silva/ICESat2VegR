#' Clips ICESat-2 ATL08 data
#'
#' @param atl08 [`ICESat2VegR::icesat2.atl08_h5-class`] object, obtained through [`ATL08_read()`]
#' for clipping
#' @param output character. Path to the output h5 file.
#' @param bbox [`numeric-class`] or [`terra::SpatExtent`] for clipping, the
#' order of the bbox is the default from NASA's ICESat-2 CMS searching:
#' (ul_lat, ul_lon, lr_lat, lr_lon).
#' @param beam [`character-class`]. The vector of beams to include, default
#' all c("gt1l", "gt2l", "gt3l", "gt1r", "gt2r", "gt3r")
#' @param additional_groups [`character-class`]. Other addional groups that should be included, default
#' c("METADATA", "orbit_info", "quality_assessment", "atlas_impulse_response", "ancillary_data")
#'
#' @return Returns the clipped S4 object of class [`ICESat2VegR::icesat2.atl08_h5-class`]
#'
#' @description This function clips the ATl08 HDF5 file. This function
#' will only clip the beam groups within hdf5, it won't change metadata
#' or ancillary data.
#'
#' @examples
#' # ATL08 file path
#' atl08_path <- system.file("extdata",
#'   "atl08_clip.h5",
#'   package = "ICESat2VegR"
#' )
#'
#' # Reading ATL08 data (h5 file)
# atl08_h5<-ATL08_read(atl08_path=atl08_path)
#'
#' # Bounding rectangle coordinates
#' xmin <- -107.7
#' xmax <- -106.5
#' ymin <- 32.75
#' ymax <- 42.75
#'
#' # Clipping ATL08 terrain and canopy attributes by boundary box
#' atl08_seg_att_dt_clip <- ATL08_h5_clipBox(atl08_h5, outdir, xmin, xmax, ymin, ymax)
#'
#' close(atl08_h5)
#' @import hdf5r
#' @export
ATL08_h5_clipBox <- function(atl08, output, bbox,
    beam = c("gt1r", "gt2r", "gt3r", "gt1l", "gt2l", "gt3l"),
    additional_groups = c("METADATA", "orbit_info", "quality_assessment", "atlas_impulse_response", "ancillary_data")) {
  ATL08_h5_clip(atl08, output, bbox, landsegmentsMask_bbox, beam, additional_groups)
}

#' Clips ICESat-2 ATL08 data
#'
#' @param atl08 [`ICESat2VegR::icesat2.atl08_h5-class`] object, obtained through [`ATL08_read()`]
#' for clipping
#' @param output character. Path to the output h5 file.
#' @param vect [`terra::SpatVector-class`] for clipping
#' @param polygon_id [`character-class`]. The attribute name used for identifying
#' the different polygons. Default is "id"
#' @param beam [`character-class`]. The vector of beams to include, default
#' all c("gt1l", "gt2l", "gt3l", "gt1r", "gt2r", "gt3r")
#' @param additional_groups [`character-class`]. Other addional groups that should be included, default
#' c("METADATA", "orbit_info", "quality_assessment", "atlas_impulse_response", "ancillary_data")
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
#' atl08_h5 <- atl08_read(atl08_path = atl08_path)
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
#' # Clipping ATL08 photons  by boundary box extent
#' atl08_photons_dt_clip <- ATL08_h5_clipGeometry(
#'   atl08_h5,
#'   output,
#'   vect,
#'   polygon_id = "id"
#' )
#'
#' close(atl08_h5)
#' @import hdf5r
#' @export
ATL08_h5_clipGeometry <- function(atl08, output, vect, polygon_id = "id",
    beam = c("gt1r", "gt2r", "gt3r", "gt1l", "gt2l", "gt3l"),
    additional_groups = c("METADATA", "orbit_info", "quality_assessment", "atlas_impulse_response", "ancillary_data")) {
  geom <- terra::union(vect)
  ATL08_h5_clip(atl08, output, clip_obj = geom, landSegmentsMask_fn = landsegmentsMask_geom, beam, additional_groups)
}

landSegments_bbox <- function(beam, bbox) {
  latitude <- longitude <- NA
  .I <- data.table::.I

  land_segments_dt <- data.table::data.table(
    latitude = beam[["land_segments/latitude"]][],
    longitude = beam[["land_segments/longitude"]][]
  )

  landSegmentsDt <- land_segments_dt[, list(latitude, longitude, .I)][
    longitude >= bbox$xmin &
      longitude <= bbox$xmax &
      latitude >= bbox$ymin &
      latitude <= bbox$ymax
  ]

  landSegmentsDt
}

landsegmentsMask_bbox <- function(beam, bbox) {
  landSegmentsBbox <- landSegments_bbox(beam, bbox)

  landSegmentsBbox$I
}

landsegmentsMask_geom <- function(beam, geom) {
  bbox <- terra::ext(geom)
  landSegmentsBbox <- landSegments_bbox(beam, bbox)

  pts <- terra::vect(
    landSegmentsBbox,
    geom = c("longitude", "latitude"),
    crs = "epsg:4326"
  )

  res <- terra::intersect(pts, geom)
  res$I
}
