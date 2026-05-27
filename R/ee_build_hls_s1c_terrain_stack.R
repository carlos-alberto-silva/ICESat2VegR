#' Build HLS, Sentinel-1C and terrain ancillary stack in Earth Engine
#'
#' @description
#' Constructs a multi-source ancillary stack in Earth Engine composed of:
#' (1) optical features from the Harmonized Landsat and Sentinel-2 (HLS) product,
#' including several vegetation indices and spectral unmixing fractions;
#' (2) C-band SAR backscatter and radar-based indices from Sentinel-1C GRD; and
#' (3) terrain variables (elevation, slope, aspect) from a DEM.
#'
#' The user must provide a temporal window via `start_date` and `end_date`.
#' The HLS image collection is filtered to this period and mosaicked using
#' `median()`. Sentinel-1C is filtered to the same period and reduced using
#' `ee$Reducer$firstNonNull()`, following the original workflow. The final
#' stack includes neighborhood statistics and GLCM texture features.
#'
#' @param x Spatial input defining the area of interest. Can be:
#'
#'   - a character path to a vector file readable by [terra::vect()] (e.g., SHP, GPKG),
#'   - a [terra::SpatVector],
#'   - a [terra::SpatRaster] (its extent will be used),
#'   - an [sf::st_sf] or [sf::st_sfc] object,
#'   - a [terra::SpatExtent] object (e.g., from [terra::ext()]),
#'   - a numeric vector of length 4 giving an extent as
#'     `c(xmin, xmax, ymin, ymax)`.
#' @param start_date Character. Start date of the temporal filter
#'   (`YYYY-MM-DD`). **Required.**
#' @param end_date Character. End date of the temporal filter
#'   (`YYYY-MM-DD`). **Required.**
#' @param cloud_max Numeric. Maximum allowed cloud coverage percentage for HLS
#'   scenes (used in the `CLOUD_COVERAGE < cloud_max` filter). Default is 10.
#' @param buffer_m Numeric. Buffer distance in meters applied to the AOI
#'   bounding box before filtering and clipping. Default is 30.
#'
#' @return
#' A named list with the following Earth Engine objects:
#' \describe{
#'   \item{aoi}{EE geometry representing the buffered AOI bounding box.}
#'   \item{hls}{HLS image with reflectance bands, vegetation indices, and
#'     spectral unmixing fractions.}
#'   \item{s1c}{Sentinel-1C image with VV/VH backscatter and radar indices.}
#'   \item{dem}{Terrain stack including elevation, slope, and aspect.}
#'   \item{stack}{Full ancillary stack combining HLS, Sentinel-1C, and terrain,
#'     including neighborhood statistics and texture metrics.}
#' }
#'
#' @details
#' This function assumes that:
#' - Earth Engine access is provided via the `ICESat2VegR` internal
#'     `ee` handle.
#' - Dataset search and catalog ID retrieval are performed by
#'     [search_datasets()] and [get_catalog_id()] from `ICESat2VegR`.
#'
#'
#' The HLS component follows the original workflow: the collection is filtered
#' by AOI and date, filtered by cloud coverage, cloud and water masks are
#' applied, and a `median()` mosaic is computed. Reflectance bands are renamed
#' (blue, green, red, nir, swir1, swir2) and a set of vegetation indices and
#' linear spectral unmixing fractions are added.
#'
#' Sentinel-1C (COPERNICUS/S1_GRD) is filtered by AOI, date, polarization
#' (VV/VH), and instrument mode ("IW"), sorted by time, and reduced using
#' `ee$Reducer$firstNonNull()`. VV and VH are scaled, and RVI and
#' co-polarization indices are derived.
#'
#' A DEM is retrieved from the NASA catalog via [search_datasets()] and
#' [get_catalog_id()], and slope and aspect are derived following the original
#' code pattern. Neighborhood statistics (mean, min, max, stdDev) are computed
#' for the HLS and DEM components using a fixed 3x3 kernel, and GLCM texture
#' metrics are computed from scaled HLS reflectance bands.
#'
#'
#' @examples
#' \dontrun{
#' # ================================================================
#' # Example 1 - Using the function
#' # ================================================================
#' res <- ee_build_hls_s1c_terrain_stack(
#'   x          = system.file("extdata", "all_boundary.shp", package = "ICESat2VegR"),
#'   start_date = "2019-04-01",
#'   end_date   = "2019-05-31"
#' )
#'
#' full_stack <- res$stack
#'
#'
#' # ================================================================
#' # Example 2 - FULL SCRIPT (no function)
#' # Users may copy/paste and modify as needed
#' # ================================================================
#'
#' library(ICESat2VegR)
#' library(terra)
#'
#' # AOI
#' geom <- terra::vect(system.file("extdata", "all_boundary.shp", package="ICESat2VegR"))
#' bbox <- terra::ext(geom)
#'
#' aoi <- ee$Geometry$BBox(
#'   west  = bbox$xmin,
#'   south = bbox$ymin,
#'   east  = bbox$xmax,
#'   north = bbox$ymax
#' )$buffer(30)
#'
#' # ---------------------------------------------------------------
#' # HLS
#' # ---------------------------------------------------------------
#' search <- search_datasets("hls")
#' id     <- get_catalog_id(search)
#'
#' cloudMask <- 2^1 + 2^2 + 2^3
#'
#' hlsMask <- function(image) image$updateMask(!(image[["Fmask"]] & cloudMask))
#' waterMask <- function(image) image$updateMask(image[["B5"]] >= 0.2)
#'
#' hls <- ee$ImageCollection(id)$
#'   filterBounds(aoi)$
#'   filterDate("2019-04-01", "2019-05-31")$
#'   filter("CLOUD_COVERAGE < 10")$
#'   map(hlsMask)$
#'   map(waterMask)$
#'   median()$
#'   clip(aoi)
#'
#' hls <- hls[["B2","B3","B4","B5","B6","B7"]]
#' names(hls) <- c("blue","green","red","nir","swir1","swir2")
#'
#' # Vegetation indices
#' nir     <- hls[["nir"]]
#' red     <- hls[["red"]]
#' blue    <- hls[["blue"]]
#' green   <- hls[["green"]]
#' swir1   <- hls[["swir1"]]
#' swir2   <- hls[["swir2"]]
#' nir_red <- nir - red
#'
#' hls[["ndvi"]] <- (nir - red) / (nir + red)
#'
#' sigma <- 1
#' knr <- exp((nir_red)^2 / (2*sigma^2))
#' hls[["kndvi"]] <- (1 - knr) / (1 + knr)
#'
#' hls[["evi"]]  <- 2.5 * nir_red / (nir + 6*red - 7.5*blue + 1)
#' hls[["savi"]] <- 1.5 * nir_red / (nir + red + 0.5)
#'
#' p1 <- 2*nir + 1
#' hls[["msavi"]] <- p1 - sqrt(p1^2 - 8*nir_red)/2
#'
#' # Spectral unmixing
#' soil  <- c(0.14, 0.16, 0.22, 0.39, 0.45, 0.27)
#' veg   <- c(0.086,0.062,0.043,0.247,0.109,0.039)
#' water <- c(0.07,0.039,0.023,0.031,0.011,0.007)
#'
#' img <- hls[[c("blue","green","red","nir","swir1","swir2")]]
#' unmixing <- img$unmix(list(soil,veg,water))
#' names(unmixing) <- c("f_soil","f_veg","f_water")
#'
#' hls <- c(hls, unmixing)
#'
#' # More indices
#' hls[["sri"]]  <- nir / red
#' hls[["ndwi"]] <- (green - nir)/(green + nir)
#' hls[["gci"]]  <- nir/green - 1
#' hls[["wdrvi"]] <- ((0.1*nir) - red)/((0.1*nir) + red)
#' hls[["gvmi"]]  <- ((nir+0.1)-(swir1+0.02))/((nir+0.1)+(swir1+0.02))
#' hls[["cvi"]]   <- nir * (red/(green^2))
#' hls[["cmr"]]   <- swir1/swir2
#'
#'
#' # ---------------------------------------------------------------
#' # DEM (NASA)
#' # ---------------------------------------------------------------
#' dem_search <- search_datasets("nasa","dem")
#' dem_id     <- get_catalog_id(dem_search)
#'
#' elevation <- ee$Image(dem_id)
#' the_slope  <- as.integer(slope(as.integer(elevation)) * 1000)
#' the_aspect <- aspect(elevation)
#'
#' stackDem <- c(elevation, the_slope, the_aspect)$clip(aoi)
#'
#'
#' # ---------------------------------------------------------------
#' # Sentinel-1C
#' # ---------------------------------------------------------------
#' s1c <- ee$ImageCollection("COPERNICUS/S1_GRD")$
#'   filterBounds(aoi)$
#'   filterDate("2019-04-01","2019-05-31")$
#'   filter(ee$Filter$listContains("transmitterReceiverPolarisation","VV"))$
#'   filter(ee$Filter$listContains("transmitterReceiverPolarisation","VH"))$
#'   filter(ee$Filter$eq("instrumentMode","IW"))
#'
#' s1c <- s1c$sort("system:time_start", FALSE)
#' s1c <- s1c$reduce(ee$Reducer$firstNonNull())
#' s1c <- s1c[["VV_first","VH_first"]]
#' names(s1c) <- c("vv","vh")
#'
#' s1c <- as.integer(s1c * 100)
#' vv <- s1c[["vv"]]
#' vh <- s1c[["vh"]]
#'
#' s1c[["rvi"]]   <- as.integer(sqrt(vv/(vv+vh)) * (vv/vh) * 10000)
#' s1c[["copol"]] <- as.integer((vv/vh) * 10000)
#' s1c[["copol2"]] <- as.integer(((vv - vh)/(vv + vh)) * 10000)
#' s1c[["copol3"]] <- as.integer((vh/vv) * 10000)
#'
#'
#' # ---------------------------------------------------------------
#' # Final stack + neighborhood + texture
#' # ---------------------------------------------------------------
#' fullStack <- c(hls, s1c, stackDem)
#'
#' kernel <- ee$Kernel$fixed(
#'   3, 3,
#'   list(c(1,1,1), c(1,1,1), c(1,1,1)),
#'   3, 3, FALSE
#' )
#'
#' for (nm in c("mean","min","max","stdDev")) {
#'   reducer <- ee$Reducer[[nm]]()
#'
#'   fullStack <- c(
#'     fullStack,
#'     hls$reduceNeighborhood(reducer, kernel),
#'     stackDem$reduceNeighborhood(reducer, kernel)
#'   )
#' }
#'
#' img_tex <- as.integer(
#'   (hls[[c("blue","green","red","nir","swir1","swir2")]]) * 1e4
#' )
#'
#' tex <- img_tex$glcmTexture(size=3)
#'
#' fullStack <- c(fullStack, tex)
#'
#' # fullStack is the ancillary layers image
#' fullStack
#' }
#'
#' @export
ee_build_hls_s1c_terrain_stack <- function(
    x,
    start_date,
    end_date,
    cloud_max = 10,
    buffer_m  = 30
) {
  if (missing(start_date) || missing(end_date)) {
    stop("'start_date' and 'end_date' are required (YYYY-MM-DD).")
  }

  # -----------------------------
  # AOI from x (path / Spat* / sf / extent / numeric)
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

    # user passed terra::ext() directly
    bbox <- x
    geom <- terra::as.polygons(bbox)

  } else if (is.numeric(x) && length(x) == 4) {

    # numeric extent vector: c(xmin, xmax, ymin, ymax)
    vec  <- unname(x)
    bbox <- terra::ext(vec[1], vec[2], vec[3], vec[4])
    geom <- terra::as.polygons(bbox)

  } else {

    stop(
      "Argument 'x' must be a path, SpatVector, SpatRaster, sf, ",
      "SpatExtent, or numeric vector of length 4 (xmin, xmax, ymin, ymax)."
    )
  }

  # Always obtain extent from geom
  bbox <- terra::ext(geom)

  # Build EE AOI from bbox
  aoi <- ee$Geometry$BBox(
    west  = bbox$xmin,
    south = bbox$ymin,
    east  = bbox$xmax,
    north = bbox$ymax
  )

  aoi <- aoi$buffer(buffer_m)

  # -----------------------------
  # HLS collection & masks
  # -----------------------------
  search <- search_datasets("hls")
  id     <- get_catalog_id(search)
  collection <- ee$ImageCollection(id)

  cloudMask <- 2^1 + 2^2 + 2^3

  hlsMask <- function(image) {
    image$updateMask(!(image[["Fmask"]] & cloudMask))
  }

  waterMask <- function(image) {
    image$updateMask(image[["B5"]] >= 0.2)
  }

  hls <- collection$
    filterBounds(aoi)$
    filterDate(start_date, end_date)$
    filter(sprintf("CLOUD_COVERAGE < %d", cloud_max))$
    map(hlsMask)$
    map(waterMask)$
    median()$
    clip(aoi)

  # Base bands
  hls <- hls[["B2", "B3", "B4", "B5", "B6", "B7"]]
  names(hls) <- c("blue", "green", "red", "nir", "swir1", "swir2")

  # -----------------------------
  # Vegetation indices
  # -----------------------------
  # NDVI
  hls[["ndvi"]] <- (hls[["nir"]] - hls[["red"]]) / (hls[["nir"]] + hls[["red"]])

  # KNDVI
  sigma   <- 1
  nir     <- hls[["nir"]]
  red     <- hls[["red"]]
  nir_red <- nir - red
  knr     <- exp((nir_red)^2 / (2 * sigma^2))
  hls[["kndvi"]] <- (1 - knr) / (1 + knr)

  # EVI
  blue <- hls[["blue"]]
  hls[["evi"]] <- 2.5 * (nir_red) / (nir + 6 * red - 7.5 * blue + 1)

  # SAVI
  hls[["savi"]] <- 1.5 * (nir_red) / (nir + red + 0.5)

  # MSAVI2
  p1 <- 2 * nir + 1
  hls[["msavi"]] <- p1 - sqrt((p1^2) - 8 * (nir_red)) / 2

  # -----------------------------
  # Linear spectral unmixing
  # -----------------------------
  soil  <- c(0.14, 0.16, 0.22, 0.39, 0.45, 0.27)
  veg   <- c(0.086, 0.062, 0.043, 0.247, 0.109, 0.039)
  water <- c(0.07, 0.039, 0.023, 0.031, 0.011, 0.007)

  img <- hls[[c("blue", "green", "red", "nir", "swir1", "swir2")]]
  unmixing <- img$unmix(list(soil, veg, water))
  names(unmixing) <- c("f_soil", "f_veg", "f_water")

  # combine HLS + unmixing
  hls <- c(hls, unmixing)

  # Simple Ratio Index
  hls[["sri"]] <- nir / red

  # NDWI
  green <- hls[["green"]]
  hls[["ndwi"]] <- (green - nir) / (green + nir)

  # Green Chlorophyll Index
  hls[["gci"]] <- (nir / green) - 1

  # WDRVI
  hls[["wdrvi"]] <- ((0.1 * nir) - red) / ((0.1 * nir) + red)

  # Global Vegetation Moisture Index
  swir1 <- hls[["swir1"]]
  hls[["gvmi"]] <- ((nir + 0.1) - (swir1 + 0.02)) / ((nir + 0.1) + (swir1 + 0.02))

  # Chlorophyll Vegetation Index
  hls[["cvi"]] <- nir * (red / (green^2))

  # Clay minerals ratio
  swir2 <- hls[["swir2"]]
  hls[["cmr"]] <- swir1 / swir2

  # -----------------------------
  # Elevation (DEM) - original pattern
  # -----------------------------
  result     <- search_datasets("nasa", "dem")
  catalog_id <- get_catalog_id(result)

  elevation <- ee$Image(catalog_id)

  the_slope  <- as.integer(slope(as.integer(elevation)) * 1000)
  the_aspect <- aspect(elevation)

  stackDem <- c(elevation, the_slope, the_aspect)$clip(aoi)

  # -----------------------------
  # Sentinel-1C - original pattern
  # -----------------------------
  s1c <- ee$ImageCollection("COPERNICUS/S1_GRD")$
    filterBounds(aoi)$
    filterDate(start_date, end_date)$
    filter(ee$Filter$listContains("transmitterReceiverPolarisation", "VV"))$
    filter(ee$Filter$listContains("transmitterReceiverPolarisation", "VH"))$
    filter(ee$Filter$eq("instrumentMode", "IW"))

  s1c <- s1c$sort("system:time_start", FALSE)
  s1c <- s1c$reduce(ee$Reducer$firstNonNull())
  s1c <- s1c[["VV_first", "VH_first"]]
  names(s1c) <- c("vv", "vh")
  s1c <- as.integer(s1c * 100)

  vv <- s1c[["vv"]]
  vh <- s1c[["vh"]]

  # RVI
  s1c[["rvi"]] <- as.integer(sqrt(vv / (vv + vh)) * (vv / vh) * 10000)

  # COPOLs
  s1c[["copol"]]  <- as.integer((vv / vh) * 10000)
  s1c[["copol2"]] <- as.integer(((vv - vh) / (vv + vh)) * 10000)
  s1c[["copol3"]] <- as.integer((vh / vv) * 10000)

  # -----------------------------
  # Final stack: HLS + S1C + DEM
  # -----------------------------
  fullStack <- c(hls, s1c, stackDem)

  kernel <- ee$Kernel$fixed(
    3, 3,
    list(c(1, 1, 1), c(1, 1, 1), c(1, 1, 1)),
    3, 3, FALSE
  )

  reducerNames <- c("mean", "min", "max", "stdDev")

  for (reducerName in reducerNames) {
    reducer <- ee$Reducer[[reducerName]]()

    hls_nb      <- hls$reduceNeighborhood(reducer, kernel)
    stackDem_nb <- stackDem$reduceNeighborhood(reducer, kernel)

    # append neighborhood stacks
    fullStack <- c(fullStack, hls_nb, stackDem_nb)
  }

  # -----------------------------
  # Texture (GLCM) from HLS bands
  # -----------------------------
  img_tex <- as.integer(
    (hls[[c("blue", "green", "red", "nir", "swir1", "swir2")]]) * 1e4
  )

  tex <- img_tex$glcmTexture(size = 3)
  fullStack <- c(fullStack, tex)

  # -----------------------------
  # Return
  # -----------------------------
  list(
    aoi   = aoi,
    hls   = hls,
    s1c   = s1c,
    dem   = stackDem,
    stack = fullStack
  )
}
