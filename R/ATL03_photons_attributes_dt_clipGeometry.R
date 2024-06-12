#' Clip ATL03 photons by Coordinates
#'
#' @description This function clips ATL03 photon attributes within given bounding coordinates.
#'
#' @param atl03_photons_dt An ATL03 photon data table. An S4 object of class [`ICESat2VegR::icesat2.atl03_dt-class`].
#' @param sppoly Spatial Polygon. An object of class [`terra::SpatVector`],
#'   which can be loaded as an ESRI shapefile using the [terra::vect] function in the
#'   \emph{sf} package.
#' @param split_by Polygon id. If defined, GEDI data will be clipped by each polygon using
#'   the polygon id from the attribute table defined by the user.
#'
#' @return Returns an S4 object of class [`ICESat2VegR::icesat2.atl03_dt-class`]
#'   containing the ATL03 photon attributes.
#'
#' @seealso \url{https://icesat-2.gsfc.nasa.gov/sites/default/files/page_files/ICESat2_ATL03_ATBD_r006.pdf}
#'
#' @examples
#' # ATL03 file path
#' atl03_path <- system.file("extdata",
#'   "atl03_clip.h5",
#'   package = "ICESat2VegR"
#' )
#'
#' # Reading ATL03 data (h5 file)
#' atl03_h5 <- ATL03_read(atl03_path = atl03_path)
#'
#' # Extracting ATL03 photon attributes
#' atl03_photons_dt <- ATL03_photons_attributes_dt(atl03_h5 = atl03_h5)
#'
#' # Specifying the path to shapefile
#' polygon_filepath <-
#'   system.file(
#'     "extdata",
#'     "clip_geom.shp",
#'     package = "ICESat2VegR"
#'   )
#'
#' # Reading shapefile as sf object
#' sppoly <- terra::vect(polygon_filepath)
#'
#' # Clipping ATL03 photon attributes by Geometry
#' atl03_photons_dt_clip <-
#'   ATL03_photons_attributes_dt_clipGeometry(atl03_photons_dt, sppoly, split_by = "id")
#'
#' head(atl03_photons_dt_clip)
#'
#' close(atl03_h5)
#' @import hdf5r stats
#' @export
ATL03_photons_attributes_dt_clipGeometry <- function(atl03_photons_dt, sppoly, split_by = "id") {
  # Check if atl03_photons_dt is of the correct class
  if (!inherits(atl03_photons_dt, "icesat2.atl03_dt")) {
    stop("atl03_photons_dt needs to be an object of class 'icesat2.atl03_dt'")
  }

  # Extract the bounding box of the spatial polygon
  exshp <- terra::ext(sppoly)

  # Clip ATL03 photons using the bounding box of the spatial polygon
  atl03_photons_dt <- ATL03_photons_attributes_dt_clipBox(
    atl03_photons_dt,
    exshp$xmin,
    exshp$xmax,
    exshp$ymin,
    exshp$ymax
  )

  # Remove any rows with NA values
  if (any(is.na(atl03_photons_dt))) {
    atl03_photons_dt <- na.omit(atl03_photons_dt)
  }

  # Add a unique identifier to each row
  atl03_photons_dt$nid <- seq_len(nrow(atl03_photons_dt))

  # Check if any data remains after clipping
  if (nrow(atl03_photons_dt) == 0) {
    warning("The polygon does not overlap the ATL03 data")
  } else {
    # Convert the ATL03 photon coordinates to a spatial vector
    points <- terra::vect(
      data.table(
        lon_ph = atl03_photons_dt$lon_ph,
        lat_ph = atl03_photons_dt$lat_ph
      ),
      geom = c("lon_ph", "lat_ph"),
      crs = terra::crs(sppoly)
    )

    # Add a row number to each point
    points$rowNumber <- as.integer(seq_along(points))

    # Intersect the points with the spatial polygon
    pts <- terra::intersect(terra::makeValid(points), terra::makeValid(sppoly))

    # If split_by is defined, clip the data by each polygon id
    if (!is.null(split_by)) {
      if (any(names(sppoly) == split_by)) {
        newFile <- atl03_photons_dt[pts$nid, ]
        newFile$poly_id <- pts[[split_by]]
      } else {
        stop("The ", split_by, " is not included in the attribute table. Please check the names in the attribute table")
      }
    } else {
      newFile <- atl03_photons_dt[pts$nid, ]
    }

    return(newFile)
  }
}
