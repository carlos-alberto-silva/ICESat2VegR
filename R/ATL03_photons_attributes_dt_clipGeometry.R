#' Clip ATL08 Terrain and Canopy Attributes by Geometry
#'
#' @description This function clips ATL08 Terrain and Canopy Attributes within a given geometry
#'
#' @param atl08_seg_att_dt A atl08_seg_att_dt object (output of
#' [`ATL08_seg_attributes_dt()`] function). An S4 object of class
#' [`ICESat2VegR::icesat2.atl08_dt-class`]
#' @param clip_obj clip_obj. An object of class [`terra::SpatVector`],
#' which can be loaded as an ESRI shapefile using [terra::vect] function in the
#' \emph{sf} package.
#' @param split_by clip_obj id. If defined, ATL08 data will be clipped by each clip_obj using
#' the clip_obj id from table of attribute defined by the user
#'
#' @return Returns an S4 object of class [`ICESat2VegR::icesat2.atl08_dt-class`]
#' containing the clipped ATL08 Terrain and Canopy Attributes.
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
#' # Extracting ATL08-derived terrain and canopy attributes
#' atl08_seg_att_dt <- ATL08_seg_attributes_dt(atl08_h5 = atl08_h5)
#'
#' clip_obj_path <- system.file("extdata",
#'   "clip_geom.shp",
#'   package = "ICESat2VegR"
#' )
#'
#' if (require(terra)) {
#'   polygon <- terra::vect(clip_obj_path)
#'
#'   head(atl08_seg_att_dt)
#'   # Clipping ATL08 Terrain and Canopy Attributes by Geometry
#'   atl08_seg_att_dt_clip <- ATL08_seg_attributes_dt_clipGeometry(
#'    atl08_seg_att_dt,
#'    polygon,
#'    split_by = "id"
#'   )
#'
#'   hasLeaflet <- require(leaflet)
#'
#'   if (hasLeaflet) {
#'     leaflet() %>%
#'       addCircleMarkers(atl08_seg_att_dt_clip$longitude,
#'         atl08_seg_att_dt_clip$latitude,
#'         radius = 1,
#'         opacity = 1,
#'         color = "red"
#'       ) %>%
#'       addScaleBar(options = list(imperial = FALSE)) %>%
#'       addPolygons(
#'         data = polygon, weight = 1, col = "white",
#'         opacity = 1, fillOpacity = 0
#'       ) %>%
#'       addProviderTiles(providers$Esri.WorldImagery,
#'         options = providerTileOptions(minZoom = 3, maxZoom = 17)
#'       )
#'   }
#' }
#' close(atl08_h5)
#' @export
ATL08_seg_attributes_dt_clipGeometry <- function(atl08_seg_att_dt, clip_obj, split_by = NULL) {
  if (!inherits(atl08_seg_att_dt, "icesat2.atl08_dt")) {
    stop("atl08_seg_att_dt needs to be an object of class 'icesat2.atl08_dt' ")
  }
  newFile <- data.table::data.table()
  exshp <- terra::ext(clip_obj)

    atl08_seg_att_dt <- ATL08_seg_attributes_dt_clipBox(
      atl08_seg_att_dt,
      exshp
    )

  if (any(is.na(atl08_seg_att_dt))) {
    atl08_seg_att_dt <- na.omit(atl08_seg_att_dt)
  } else {
    atl08_seg_att_dt <- atl08_seg_att_dt
  }

  atl08_seg_att_dt$nid <- 1:nrow(atl08_seg_att_dt)

  if (nrow(atl08_seg_att_dt) == 0) {
    message("The clip_obj does not overlap the ATL08 data")
  } else {
    points <- terra::vect(
      as.data.frame(atl08_seg_att_dt),
      geom = c("longitude", "latitude"),
      crs = terra::crs(clip_obj)
    )

    points$rowNumber <- as.integer(seq_along(points))
    pts <- terra::intersect(terra::makeValid(points), terra::makeValid(clip_obj))

    if (!is.null(split_by)) {
      if (any(names(clip_obj) == split_by)) {
        newFile <- atl08_seg_att_dt[pts$nid, ]
        newFile$poly_id <- pts[[split_by]]
      } else {
        stop(paste("The", split_by, "is not included in the attribute table.
                       Please check the names in the attribute table"))
      }
    } else {
      newFile <- atl08_seg_att_dt[pts$nid, ]
    }

    prepend_class(newFile, "icesat2.atl08_dt")

  }
  return(newFile)
}
