#' Build HLS, Sentinel-1C and terrain ancillary stack in Earth Engine
#'
#' @description
#' Constructs a multi-source ancillary stack in Google Earth Engine composed of:
#' \enumerate{
#'   \item Optical features from the Harmonized Landsat and Sentinel-2 (HLS)
#'     product, including spectral bands, vegetation indices, and spectral
#'     unmixing fractions.
#'   \item C-band SAR backscatter and radar-based indices from Sentinel-1C GRD.
#'   \item Terrain variables (elevation, slope, aspect) from a NASA DEM.
#'   \item Neighborhood statistics (mean, min, max, stdDev) computed over a
#'     3x3 kernel for HLS and terrain bands.
#'   \item GLCM texture metrics computed from scaled HLS reflectance bands.
#' }
#'
#' The final stack contains 244 bands covering spectral, radar, terrain,
#' neighborhood, and texture features.
#'
#' \strong{Authentication:} This function requires an active Google Earth
#' Engine session. Authenticate using \code{ICESat2VegR::ee_initialize()}
#' before calling this function.
#'
#' @param x Spatial input defining the area of interest. Can be:
#'   \itemize{
#'     \item a character path to a vector file readable by
#'       \code{\link[terra]{vect}} (e.g., SHP, GPKG),
#'     \item a \code{\link[terra]{SpatVector}},
#'     \item a \code{\link[terra]{SpatRaster}} (its extent will be used),
#'     \item an \code{sf} or \code{sfc} object,
#'     \item a \code{\link[terra]{SpatExtent}} object,
#'     \item a numeric vector of length 4:
#'       \code{c(xmin, xmax, ymin, ymax)}.
#'   }
#' @param start_date character. Start date of the temporal filter
#'   in \code{YYYY-MM-DD} format. \strong{Required.}
#' @param end_date character. End date of the temporal filter
#'   in \code{YYYY-MM-DD} format. \strong{Required.}
#' @param cloud_max numeric. Maximum allowed cloud coverage percentage
#'   for HLS scenes. Default is \code{10}.
#' @param buffer_m numeric. Buffer distance in meters applied to the AOI
#'   bounding box before filtering and clipping. Default is \code{30}.
#' @param hls_collection_id character vector. HLS collection IDs to try
#'   in order until one has data for the AOI and date range. Default tries
#'   Sentinel-2 HLS first, then Landsat HLS as fallback:
#'   \code{c("NASA/HLS/HLSS30/v002", "NASA/HLS/HLSL30/v002")}.
#'
#' @return Returns an \code{ee$Image} object containing 244 bands:
#'   \itemize{
#'     \item HLS reflectance bands and vegetation indices
#'       (blue, green, red, nir, swir1, swir2, ndvi, kndvi, evi,
#'       savi, msavi, f_soil, f_veg, f_water, sri, ndwi, gci,
#'       wdrvi, gvmi, cvi, cmr).
#'     \item Sentinel-1C SAR bands and radar indices
#'       (vv, vh, rvi, copol, copol2, copol3).
#'     \item Terrain bands (elevation, slope, aspect).
#'     \item Neighborhood statistics (mean, min, max, stdDev)
#'       for HLS and terrain bands.
#'     \item GLCM texture metrics for HLS reflectance bands.
#'   }
#'
#' @details
#' The HLS collection is selected automatically from \code{hls_collection_id}
#' by trying each ID in order and using the first one that contains scenes
#' for the specified AOI and date range. By default, Sentinel-2 HLS
#' (\code{NASA/HLS/HLSS30/v002}) is tried first, followed by Landsat HLS
#' (\code{NASA/HLS/HLSL30/v002}).
#'
#' Cloud masks are applied using the \code{Fmask} band before computing
#' the median mosaic. Sentinel-1C is filtered by VV/VH polarization and
#' IW instrument mode, sorted by acquisition time, and reduced using
#' \code{ee$Reducer$firstNonNull()}.
#'
#' Neighborhood statistics (mean, min, max, stdDev) are computed over a
#' fixed 3x3 kernel for HLS and terrain bands. GLCM texture metrics are
#' computed from scaled (x10000) HLS reflectance bands using a window
#' size of 3.
#'
#' @examples
#' \dontrun{
#'   # Requires Google Earth Engine authentication
#'   Sys.setenv(EE_PROJECT = "your-ee-project-id")
#'   ICESat2VegR::ee_initialize()
#'
#'   library(sf)
#'   library(terra)
#'
#'   # AOI in Ocala National Forest, Florida
#'   aoi <- st_as_sfc(st_bbox(c(
#'     xmin = -81.8, xmax = -81.6,
#'     ymin =  29.1, ymax =  29.3
#'   ), crs = 4326))
#'
#'   # Build the full ancillary stack (244 bands)
#'   stack <- ee_build_hls_s1c_terrain_stack(
#'     x          = aoi,
#'     start_date = "2024-03-01",
#'     end_date   = "2024-06-01",
#'     cloud_max  = 10
#'   )
#'
#'   # Check available bands
#'   names(stack)
#'
#'   # Compute AOI centroid for map view
#'   centroid_lon <- mean(terra::ext(terra::vect(aoi))[c(1, 2)])
#'   centroid_lat <- mean(terra::ext(terra::vect(aoi))[c(3, 4)])
#'
#'   # Visualize true color RGB
#'   leaflet::leaflet() |>
#'     leaflet::addProviderTiles("Esri.WorldImagery") |>
#'     ICESat2VegR::addEEImage(
#'       stack$select(c("red", "green", "blue")),
#'       bands = c("red", "green", "blue"),
#'       group = "True Color",
#'       min   = 0,
#'       max   = 0.3
#'     ) |>
#'     leaflet::addControl(
#'       html     = "<b>True Color</b><br>red / green / blue",
#'       position = "bottomleft"
#'     ) |>
#'     leaflet::setView(lng = centroid_lon, lat = centroid_lat, zoom = 11) |>
#'     leaflet::addLayersControl(
#'       overlayGroups = "True Color",
#'       options = leaflet::layersControlOptions(collapsed = FALSE)
#'     )
#'
#'   # Visualize false color (nir, red, green)
#'   leaflet::leaflet() |>
#'     leaflet::addProviderTiles("Esri.WorldImagery") |>
#'     ICESat2VegR::addEEImage(
#'       stack$select(c("nir", "red", "green")),
#'       bands = c("nir", "red", "green"),
#'       group = "False Color",
#'       min   = 0,
#'       max   = 0.3
#'     ) |>
#'     leaflet::addControl(
#'       html     = "<b>False Color</b><br>nir / red / green",
#'       position = "bottomleft"
#'     ) |>
#'     leaflet::setView(lng = centroid_lon, lat = centroid_lat, zoom = 11) |>
#'     leaflet::addLayersControl(
#'       overlayGroups = "False Color",
#'       options = leaflet::layersControlOptions(collapsed = FALSE)
#'     )
#'
#'   # Visualize NDVI
#'   leaflet::leaflet() |>
#'     leaflet::addProviderTiles("Esri.WorldImagery") |>
#'     ICESat2VegR::addEEImage(
#'       stack$select("ndvi"),
#'       bands = "ndvi",
#'       group = "NDVI",
#'       min   = 0,
#'       max   = 0.8
#'     ) |>
#'     leaflet::addControl(
#'       html     = "<b>NDVI</b>",
#'       position = "bottomleft"
#'     ) |>
#'     leaflet::setView(lng = centroid_lon, lat = centroid_lat, zoom = 11) |>
#'     leaflet::addLayersControl(
#'       overlayGroups = "NDVI",
#'       options = leaflet::layersControlOptions(collapsed = FALSE)
#'     )
#'
#'   # Force Landsat HLS only
#'   stack_landsat <- ee_build_hls_s1c_terrain_stack(
#'     x                 = aoi,
#'     start_date        = "2024-03-01",
#'     end_date          = "2024-06-01",
#'     hls_collection_id = "NASA/HLS/HLSL30/v002"
#'   )
#'   names(stack_landsat)
#' }
#'
#' @export
ee_build_hls_s1c_terrain_stack <- function(
    x,
    start_date,
    end_date,
    cloud_max = 10,
    buffer_m  = 30,
    hls_collection_id = c("NASA/HLS/HLSS30/v002", "NASA/HLS/HLSL30/v002")
) {
  if (missing(start_date) || missing(end_date)) {
    stop("'start_date' and 'end_date' are required (YYYY-MM-DD).")
  }

  # -----------------------------
  # AOI from x
  # -----------------------------
  if (is.character(x)) {
    geom <- terra::vect(x)
  } else if (inherits(x, "SpatVector")) {
    geom <- x
  } else if (inherits(x, "SpatRaster")) {
    geom <- terra::as.polygons(terra::ext(x))
  } else if (inherits(x, "sf") || inherits(x, "sfc")) {
    geom <- terra::vect(x)
  } else if (inherits(x, "SpatExtent")) {
    geom <- terra::as.polygons(x)
  } else if (is.numeric(x) && length(x) == 4) {
    vec  <- unname(x)
    geom <- terra::as.polygons(terra::ext(vec[1], vec[2], vec[3], vec[4]))
  } else {
    stop(
      "Argument 'x' must be a path, SpatVector, SpatRaster, sf, ",
      "SpatExtent, or numeric vector of length 4 (xmin, xmax, ymin, ymax)."
    )
  }

  bbox <- terra::ext(geom)
  aoi  <- ee$Geometry$BBox(
    west  = bbox$xmin,
    south = bbox$ymin,
    east  = bbox$xmax,
    north = bbox$ymax
  )$buffer(buffer_m)

  # -----------------------------
  # HLS collection
  # -----------------------------
  hls_collection <- NULL
  for (hls_id in hls_collection_id) {
    candidate <- ee$ImageCollection(hls_id)$
      filterBounds(aoi)$
      filterDate(start_date, end_date)$
      filter(sprintf("CLOUD_COVERAGE < %d", cloud_max))
    count <- candidate$size()$getInfo()
    if (count > 0) {
      hls_collection <- candidate
      message(sprintf("Using HLS collection: %s (%d scenes found)", hls_id, count))
      break
    }
  }

  if (is.null(hls_collection)) {
    stop(
      "No HLS scenes found for the specified AOI and date range. ",
      "Try expanding the date range or increasing cloud_max."
    )
  }

  cloudMask <- 2^1 + 2^2 + 2^3
  hlsMask   <- function(image) image$updateMask(!(image[["Fmask"]] & cloudMask))
  #waterMask <- function(image) image$updateMask(image[["B5"]] >= 0.2)

  hls <- hls_collection$
    map(hlsMask)$
    #map(waterMask)$
    median()$
    clip(aoi)

  hls <- hls[["B2", "B3", "B4", "B5", "B6", "B7"]]
  names(hls) <- c("blue", "green", "red", "nir", "swir1", "swir2")

  # -----------------------------
  # Vegetation indices
  # -----------------------------
  hls[["ndvi"]] <- (hls[["nir"]] - hls[["red"]]) / (hls[["nir"]] + hls[["red"]])

  sigma   <- 1
  nir     <- hls[["nir"]]
  red     <- hls[["red"]]
  nir_red <- nir - red
  knr     <- exp((nir_red)^2 / (2 * sigma^2))
  hls[["kndvi"]] <- (1 - knr) / (1 + knr)

  blue <- hls[["blue"]]
  hls[["evi"]]   <- 2.5 * (nir_red) / (nir + 6 * red - 7.5 * blue + 1)
  hls[["savi"]]  <- 1.5 * (nir_red) / (nir + red + 0.5)

  p1 <- 2 * nir + 1
  hls[["msavi"]] <- p1 - sqrt((p1^2) - 8 * (nir_red)) / 2

  # -----------------------------
  # Spectral unmixing
  # -----------------------------
  soil  <- c(0.14, 0.16, 0.22, 0.39, 0.45, 0.27)
  veg   <- c(0.086, 0.062, 0.043, 0.247, 0.109, 0.039)
  water <- c(0.07, 0.039, 0.023, 0.031, 0.011, 0.007)

  img      <- hls[[c("blue", "green", "red", "nir", "swir1", "swir2")]]
  unmixing <- img$unmix(list(soil, veg, water))
  names(unmixing) <- c("f_soil", "f_veg", "f_water")
  hls <- c(hls, unmixing)

  # More indices
  hls[["sri"]]   <- nir / red
  green <- hls[["green"]]
  hls[["ndwi"]]  <- (green - nir) / (green + nir)
  hls[["gci"]]   <- (nir / green) - 1
  hls[["wdrvi"]] <- ((0.1 * nir) - red) / ((0.1 * nir) + red)
  swir1 <- hls[["swir1"]]
  hls[["gvmi"]]  <- ((nir + 0.1) - (swir1 + 0.02)) / ((nir + 0.1) + (swir1 + 0.02))
  hls[["cvi"]]   <- nir * (red / (green^2))
  swir2 <- hls[["swir2"]]
  hls[["cmr"]]   <- swir1 / swir2

  # -----------------------------
  # DEM
  # -----------------------------
  dem_result <- search_datasets("nasa", "dem")
  catalog_id <- get_catalog_id(dem_result)
  elevation  <- ee$Image(catalog_id)
  the_slope  <- as.integer(slope(as.integer(elevation)) * 1000)
  the_aspect <- aspect(elevation)
  stackDem   <- c(elevation, the_slope, the_aspect)$clip(aoi)
  # -----------------------------
  # Sentinel-1C
  # -----------------------------
  s1c <- ee$ImageCollection("COPERNICUS/S1_GRD")$
    filterBounds(aoi)$
    filterDate(start_date, end_date)$
    filter(ee$Filter$listContains("transmitterReceiverPolarisation", "VV"))$
    filter(ee$Filter$listContains("transmitterReceiverPolarisation", "VH"))$
    filter(ee$Filter$eq("instrumentMode", "IW"))

  s1c <- s1c$sort("system:time_start", FALSE)$
    reduce(ee$Reducer$firstNonNull())
  s1c <- s1c[["VV_first", "VH_first"]]
  names(s1c) <- c("vv", "vh")
  s1c <- as.integer(s1c * 100)

  vv <- s1c[["vv"]]
  vh <- s1c[["vh"]]
  s1c[["rvi"]]    <- as.integer(sqrt(vv / (vv + vh)) * (vv / vh) * 10000)
  s1c[["copol"]]  <- as.integer((vv / vh) * 10000)
  s1c[["copol2"]] <- as.integer(((vv - vh) / (vv + vh)) * 10000)
  s1c[["copol3"]] <- as.integer((vh / vv) * 10000)

  # -----------------------------
  # Full stack
  # -----------------------------
  fullStack <- c(hls, s1c, stackDem)

  kernel <- ee$Kernel$fixed(
    3, 3,
    list(c(1, 1, 1), c(1, 1, 1), c(1, 1, 1)),
    3, 3, FALSE
  )

  for (reducerName in c("mean", "min", "max", "stdDev")) {
    reducer     <- ee$Reducer[[reducerName]]()
    hls_nb      <- hls$reduceNeighborhood(reducer, kernel)
    stackDem_nb <- stackDem$reduceNeighborhood(reducer, kernel)
    fullStack   <- c(fullStack, hls_nb, stackDem_nb)
  }

  img_tex   <- as.integer((hls[[c("blue", "green", "red", "nir", "swir1", "swir2")]]) * 1e4)
  tex       <- img_tex$glcmTexture(size = 3)
  fullStack <- c(fullStack, tex)

  # -----------------------------
  # Return fullStack only
  # -----------------------------
  invisible(fullStack)
}
