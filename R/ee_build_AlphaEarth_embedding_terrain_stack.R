#' Stack Alpha Earth embedding and terrain ancillary layers
#'
#' @description
#' Creates a Google Earth Engine image stack that combines:
#'
#'   - Annual satellite embeddings from
#'     \code{GOOGLE/SATELLITE_EMBEDDING/V1/ANNUAL} (AlphaEarth).
#'   - Terrain derivatives (elevation, slope, aspect) from USGS 3DEP or
#'     NASADEM as fallbacks.
#'   - Optional longitude/latitude bands from \code{ee$Image$pixelLonLat()}.
#'
#' The embedding collection is filtered to the specified year range and AOI,
#' reduced using a pixel-wise median, and clipped to the AOI. Terrain data
#' are reprojected to a target scale (default \code{30} m).
#'
#' \strong{Geographic coverage limitation:} The primary terrain sources
#' (USGS 3DEP 1 m and 10 m) cover only the contiguous United States,
#' Alaska, and Hawaii. For areas outside the US, the function falls back
#' to NASADEM, which provides global coverage but at coarser resolution
#' (~30 m). For best terrain accuracy, this function is recommended for
#' use within the United States.
#'
#' \strong{Authentication:} This function requires an active Google Earth
#' Engine session. Authenticate using \code{ICESat2VegR::ee_initialize()}
#' before calling this function.
#'
#' @param geom AOI geometry. Can be:
#'   \itemize{
#'     \item an \code{sf} or \code{sfc} object,
#'     \item a \code{\link[terra]{SpatVector}} object,
#'     \item or an Earth Engine geometry (Python object) already in
#'       geographic coordinates.
#'   }
#'   The geometry is internally converted to an Earth Engine geometry
#'   via an internal helper.
#' @param start_year,end_year Integer (or coercible to integer). Inclusive
#'   year range used to filter the AlphaEarth annual embedding collection.
#' @param mask_outside Logical. If \code{TRUE} (default), pixels outside
#'   the AOI are masked out using a binary mask image clipped to
#'   \code{geom}.
#' @param terrain_scale Numeric. Target scale in meters used to reproject
#'   the terrain data. Default is \code{30}.
#' @param multiply_slope_aspect_by10 Logical. If \code{TRUE} (default),
#'   slope and aspect are multiplied by 10 to preserve conventions used
#'   elsewhere in the package and avoid small floating-point values.
#' @param add_lonlat Logical. If \code{TRUE} (default), longitude and
#'   latitude bands named \code{lon} and \code{lat} are added, derived
#'   from \code{ee$Image$pixelLonLat()}.
#'
#' @return
#' A Python \code{ee$Image} object containing:
#'   \itemize{
#'     \item The median annual AlphaEarth embedding bands.
#'     \item \code{elevation}, \code{slope}, and \code{aspect} bands.
#'     \item Optional \code{lon} and \code{lat} bands.
#'   }
#'
#' @details
#' The terrain source is chosen in the following order:
#' \enumerate{
#'   \item USGS 3DEP 1 m (\code{USGS/3DEP/1m}), if available for the AOI.
#'   \item USGS 3DEP 10 m (\code{USGS/3DEP/10m}).
#'   \item NASADEM (\code{NASA/NASADEM_HGT/001}) as a global fallback.
#' }
#'
#' All terrain layers are clipped to the AOI and reprojected to
#' \code{EPSG:4326} with the requested \code{terrain_scale}.
#'
#' @examples
#' \dontrun{
#'   # Requires Google Earth Engine authentication
#'   # Note: USGS 3DEP terrain data covers the United States only.
#'   # Outside the US, NASADEM is used as a global fallback.
#'   Sys.setenv(EE_PROJECT = "your-ee-project-id")
#'   ICESat2VegR::ee_initialize()
#'
#'   library(sf)
#'   aoi <- st_as_sfc(st_bbox(c(
#'     xmin = -82.4, xmax = -82.2,
#'     ymin =  29.6, ymax =  29.8
#'   ), crs = 4326))
#'
#'   ee_geom <- ICESat2VegR:::.as_ee_geom(aoi)
#'
#'   # Build the embedding + terrain stack for 2025
#'   predictors <- ee_build_AlphaEarth_embedding_terrain_stack(
#'     geom       = ee_geom,
#'     start_year = 2025,
#'     end_year   = 2025
#'   )
#'
#'   # Visualize RGB embedding bands using leaflet
#'   predictors_rgb <- predictors$select(c("A00", "A20", "A40"))
#'
#'   centroid_lon <- mean(terra::ext(terra::vect(aoi))[c(1, 2)])
#'   centroid_lat <- mean(terra::ext(terra::vect(aoi))[c(3, 4)])
#'
#'   leaflet::leaflet() |>
#'     ICESat2VegR::addEEImage(
#'       predictors_rgb,
#'       bands = c("A00", "A20", "A40"),
#'       group = "RGB Embedding",
#'       min   = -0.2,
#'       max   =  0.2
#'     ) |>
#'     leaflet::setView(
#'       lng  = centroid_lon,
#'       lat  = centroid_lat,
#'       zoom = 15
#'     ) |>
#'     leaflet::addLayersControl(
#'       overlayGroups = "RGB Embedding",
#'       options = leaflet::layersControlOptions(collapsed = FALSE)
#'     )
#' }
#'
#' @export
ee_build_AlphaEarth_embedding_terrain_stack <- function(
    geom,
    start_year,
    end_year,
    mask_outside              = TRUE,
    terrain_scale             = 30,
    multiply_slope_aspect_by10 = TRUE,
    add_lonlat                = TRUE
) {
  ee <- reticulate::import("ee", delay_load = FALSE)         # Import Earth Engine
  .ee_ping(ee)                                 # Ping EE session (internal helper)

  start_year <- as.integer(start_year)
  end_year   <- as.integer(end_year)

  if (!is.finite(start_year) || !is.finite(end_year) || end_year < start_year) {
    stop("Invalid year range: 'start_year' and 'end_year' must be finite and end_year >= start_year.")
  }

  # Convert R geometry to EE geometry (internal helper handles sf/terra/EE)
  ee_geom <- .as_ee_geom(geom)

  # AlphaEarth annual embeddings (GOOGLE/SATELLITE_EMBEDDING/V1/ANNUAL)
  emb_ic <- ee$ImageCollection("GOOGLE/SATELLITE_EMBEDDING/V1/ANNUAL")$
    filterDate(
      sprintf("%04d-01-01", start_year),
      sprintf("%04d-12-31", end_year)
    )$
    filterBounds(ee_geom)

  embedding <- emb_ic$median()$clip(ee_geom)

  # Terrain: try USGS 3DEP 1 m, then 10 m, then NASADEM
  elevation <- tryCatch({
    ee$ImageCollection("USGS/3DEP/1m")$mosaic()$select("elevation")$clip(ee_geom)
  }, error = function(e) NULL)

  if (is.null(elevation)) {
    elevation <- tryCatch({
      ee$ImageCollection("USGS/3DEP/10m")$mosaic()$select("elevation")$clip(ee_geom)
    }, error = function(e) NULL)
  }

  if (is.null(elevation)) {
    elevation <- ee$Image("NASA/NASADEM_HGT/001")$select("elevation")$clip(ee_geom)
  }

  elev_proj <- elevation$reproject(crs = "EPSG:4326", scale = terrain_scale)

  slope  <- ee$Terrain$slope(elev_proj)$rename("slope")
  aspect <- ee$Terrain$aspect(elev_proj)$rename("aspect")

  if (isTRUE(multiply_slope_aspect_by10)) {
    slope  <- slope$multiply(10)
    aspect <- aspect$multiply(10)
  }

  slope  <- slope$clip(ee_geom)
  aspect <- aspect$clip(ee_geom)

  if (isTRUE(add_lonlat)) {
    lonlat <- ee$Image$pixelLonLat()$clip(ee_geom)
    lon    <- lonlat$select("longitude")$rename("lon")
    lat    <- lonlat$select("latitude")$rename("lat")
  }

  stack <- embedding$
    addBands(elevation$rename("elevation"))$
    addBands(slope)$
    addBands(aspect)

  if (isTRUE(add_lonlat)) {
    stack <- stack$addBands(lon)$addBands(lat)
  }

  if (isTRUE(mask_outside)) {
    mask  <- ee$Image$constant(1)$clip(ee_geom)
    stack <- stack$updateMask(mask)
  }

  stack
}
