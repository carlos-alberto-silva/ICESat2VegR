#' Clip ATL03 photons by Coordinates
#'
#' @description This function clips ATL03 photons attributes within a given bounding coordinates
#'
#' @usage ATL03_photons_attributes_dt_clipGeometry(atl03_photons_dt, xmin, xmax, ymin, ymax)
#'
#' @param atl03_photons_dt A atl03_photons_dt object (output of [atl03_photons_attributes_dt()] function).
#' An S4 object of class [ICESat2VegR::icesat2.atl03_dt]
#' @param sppoly Spatial Polygon. An object of class [`terra::SpatVector`],
#' which can be loaded as an ESRI shapefile using [terra::vect] function in the
#' \emph{sf} package.
#' @param split_by Polygon id. If defined, GEDI data will be clipped by each polygon using
#' the polygon id from table of attribute defined by the user
#'
#' @return Returns an S4 object of class [ICESat2VegR::icesat2.atl03_dt]
#' containing the ATL03 photons attributes.
#'
#'@seealso \url{https://icesat-2.gsfc.nasa.gov/sites/default/files/page_files/ICESat2_atl03_ATBD_r006.pdf}
#'
#' @examples
#'# Specifying the path to ATL03 file (zip file)
#'outdir = tempdir()
#'atl03_zip <- system.file("extdata",
#'                   "atl03_20220401221822_01501506_005_01.zip",
#'                   package="ICESat2VegR")
#'
#'# Unzipping ATL03 file
#'atl03_path <- unzip(atl03_zip,exdir = outdir)
#'
#'# Reading ATL03 data (h5 file)
#atl03_h5<-atl03_read(atl03_path=atl03_path)
#'
#'# Extracting ATL03 photons attributes
#'atl03_photons_dt<-ATL03_photons_attributes_dt(atl03_h5=atl03_h5)
#'
#' # Specifying the path to shapefile
#' polygon_filepath <- system.file("extdata", "polygon.shp", package = "ICESat2VegR")
#'
#' # Reading shapefile as sf object
#' sppoly <- terra::vect(polygon_filepath)
#'
#' # Clipping ATL03 photons attributes by Geometry
#' atl03_photons_dt_clip <- ATL03_photons_attributes_dt_clipGeometry(atl03_photons_dt, sppoly, split_by = "FID")
#' head(atl03_photons_dt_clip)
#'
#'close(atl03_h5)
#'@import hdf5r stats
#'@export
ATL03_photons_attributes_dt_clipGeometry <- function(atl03_photons_dt, sppoly, split_by = "id") {

  if (!inherits(atl03_photons_dt, "icesat2.atl03_dt")){
    stop("atl03_photons_dt needs to be an object of class 'icesat2.at03_dt' ")
  }

  exshp <- terra::ext(sppoly)

  atl03_photons_dt <- ATL03_photons_attributes_dt_clipBox(
    atl03_photons_dt,
    exshp$xmin,
    exshp$xmax,
    exshp$ymin,
    exshp$ymax
  )

  if (any(is.na(atl03_photons_dt))) {
    atl03_photons_dt<-na.omit(atl03_photons_dt)
  }

  atl03_photons_dt$nid<-1:nrow(atl03_photons_dt)

  if (nrow(atl03_photons_dt) == 0) {
    print("The polygon does not overlap the ATL08 data")
  } else {
    points <- terra::vect(
      data.table(lon_ph=atl03_photons_dt$lon_ph,
                 lat_ph=atl03_photons_dt$lat_ph),
      geom = c("lon_ph", "lat_ph"),
      crs = terra::crs(sppoly)
    )

    points$rowNumber <- as.integer(seq_along(points))
    pts <- terra::intersect(terra::makeValid(points), terra::makeValid(sppoly))

    if (!is.null(split_by)) {
      if (any(names(sppoly) == split_by)) {
        newFile <- atl03_photons_dt[pts$nid, ]
        newFile$poly_id<-pts[[split_by]]
      } else {
        stop(paste("The", split_by, "is not included in the attribute table.
                       Please check the names in the attribute table"))
      }
    } else {

      newFile <- atl03_photons_dt[pts$nid, ]
      #newFile <- atl03_photons_dt[mask, ]
    }

    return(newFile)
  }
}

