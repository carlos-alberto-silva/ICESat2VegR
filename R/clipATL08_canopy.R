#' Clip ATL08 Canopy Height Metrics by Coordinates
#'
#' @description This function clips ATL08 Canopy Height Metrics within a given bounding coordinates
#'
#' @usage clipATL08_canopy_box(canopy_metrics, xmin, xmax, ymin, ymax)
#'
#' @param canopy_metrics A canopy_metrics object (output of [ATL08_canopy()] function).
#' An S4 object of class [data.table::data.table]
#' @param xmin Numeric. West longitude (x) coordinate of bounding rectangle, in decimal degrees.
#' @param xmax Numeric. East longitude (x) coordinate of bounding rectangle, in decimal degrees.
#' @param ymin Numeric. South latitude (y) coordinate of bounding rectangle, in decimal degrees.
#' @param ymax Numeric. North latitude (y) coordinate of bounding rectangle, in decimal degrees.
#'
#' @return Returns an S4 object of class [data.table::data.table]
#' containing the ATL08 Canopy Height Metrics.
#'
#'@seealso \url{https://icesat-2.gsfc.nasa.gov/sites/default/files/page_files/ICESat2_ATL08_ATBD_r006.pdf}
#'
#' @examples
#'# Specifying the path to ATL08 file (zip file)
#'outdir = tempdir()
#'atl08_zip <- system.file("extdata",
#'                   "ATL08_20220401221822_01501506_005_01.zip",
#'                   package="rICESat2Veg")
#'
#'# Unzipping ATL08 file
#'atl08_path <- unzip(atl08_zip,exdir = outdir)
#'
#'# Reading ATL08 data (h5 file)
#atl08_h5<-readATL08(ATL08path=atl08_path)
#'
#'# Extracting ATL08-derived Canopy Metrics
#'canopy_metrics<-ATL08_canopy(atl08_h5=atl08_h5)
#'
#' # Bounding rectangle coordinates
#' xmin <- -107.7
#' xmax <- -106.5
#' ymin <- 32.75
#' ymax <- 42.75
#'
#' # Clipping ATL08-derived canopy metrics by boundary box extent
#'canopy_metrics_clip <- clipATL08_canopy_box(canopy_metrics, xmin, xmax, ymin, ymax)
#'
#'close(level2a)
#'@import hdf5r stats
#'@export
clipATL08_canopy_box <- function(canopy_metrics, xmin, xmax, ymin, ymax) {
  # xmin ymin xmax ymax
  mask <-
    canopy_metrics$longitude >= xmin &
    canopy_metrics$longitude <= xmax &
    canopy_metrics$latitude >= ymin &
    canopy_metrics$latitude <= ymax

  mask[!stats::complete.cases(mask)] <- FALSE
  mask <- (seq_along(canopy_metrics$latitude))[mask]
  newFile <- canopy_metrics[mask, ]
  # newFile<- new("gedi.level1b.dt", dt = level1bdt[mask,])
  return(newFile)
}

#' Clip ATL08 Canopy Height Metrics by Geometry
#'
#' @description This function clips ATL08 Canopy Height Metrics within a given geometry
#'
#' @usage clipATL08_canopy_geometry(canopy_metrics, xmin, xmax, ymin, ymax)
#'
#' @param canopy_metrics A canopy_metrics object (output of [ATL08_canopy()] function).
#' An S4 object of class "gedi.level2a".
#' @param polygon Polygon. An object of class [`sf::sf`],
#' which can be loaded as an ESRI shapefile using [sf::st_read] function in the
#' \emph{sf} package.
#' @param split_by Polygon id. If defined, GEDI data will be clipped by each polygon using
#' the polygon id from table of attribute defined by the user
#'
#' @return Returns an S4 object of class [data.table::data.table]
#' containing the clipped ATL08 Canopy Height Metrics.
#'
#' @examples
#' # Specifying the path to ATL08 file (zip file)
#'outdir = tempdir()
#'atl08_zip <- system.file("extdata",
#'                   "ATL08_20220401221822_01501506_005_01.zip",
#'                   package="rICESat2Veg")
#'
#' # Unzipping ATL08 file
#'atl08_path <- unzip(atl08_zip,exdir = outdir)
#'
#' # Reading ATL08 data (h5 file)
#atl08_h5<-readATL08(ATL08path=atl08_path)
#'
#'# Extracting ATL08-derived Canopy Metrics
#'canopy_metrics<-ATL08_canopy(atl08_h5=atl08_h5)
#'
#' # Specifying the path to shapefile
#' polygon_filepath <- system.file("extdata", "polygon.shp", package = "rICESat2Veg")
#'
#' # Reading shapefile as sf object
#' library(sf)
#' polygon <- sf::st_read(polygon_filepath)
#'
#' # Clipping ATL08 Canopy Height Metrics by Geometry
#' canopy_metrics_clip <- clipATL08_canopy_geometry(canopy_metrics, polygon, split_by = "FID")
#'
#' hasLeaflet <- require(leaflet)
#'
#' if (hasLeaflet) {
#'   leaflet() %>%
#'     addCircleMarkers(canopy_metrics_clip$longitude,
#'       canopy_metrics_clip$latitude,
#'       radius = 1,
#'       opacity = 1,
#'       color = "red"
#'     ) %>%
#'     addScaleBar(options = list(imperial = FALSE)) %>%
#'     addPolygons(
#'       data = polygon, weight = 1, col = "white",
#'       opacity = 1, fillOpacity = 0
#'     ) %>%
#'     addProviderTiles(providers$Esri.WorldImagery)
#' }
#' close(atl08_h5)
#' @export
clipATL08_canopy_geometry <- function(canopy_metrics, polygon, split_by = "id") {

  canopy_metrics$nid<-1:nrow(canopy_metrics)
  exshp <- sf::st_bbox(polygon)
  canopy_metricsdt <- clipATL08_canopy_box(
    canopy_metrics,
    xmin = exshp$xmin,
    xmax = exshp$xmax,
    ymin = exshp$ymin,
    ymax = exshp$ymax
  )
  if (nrow(canopy_metricsdt) == 0) {
    print("The polygon does not overlap the ATL08 data")
  } else {

    points <- sf::st_as_sf(
      canopy_metricsdt,
      coords = c("longitude", "latitude"),
      crs = sf::st_crs(polygon)
    )

    #names(points) <- gsub("^(?!geometry)", "x_\\1", names(points), perl = TRUE)
    pts <- sf::st_intersection(sf::st_make_valid(points), sf::st_make_valid(polygon))

    if (!is.null(split_by)) {
      if (any(names(polygon) == split_by)) {
        newFile <- canopy_metricsdt[pts$nid, ]
        newFile$poly_id<-pts[[split_by]]
      } else {
        stop(paste("The", split_by, "is not included in the attribute table.
                       Please check the names in the attribute table"))
      }
    } else {
      newFile <- canopy_metricsdt[pts$nid, ]
    }
    # newFile<- new("gedi.level1b.dt", dt = level2adt2@dt[mask,])
    return(newFile)
  }
}
