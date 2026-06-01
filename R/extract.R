#' Extract ancillary stack values at point locations from Google Earth Engine
#'
#' @description
#' Extracts pixel values from a Google Earth Engine image stack at a set of
#' point locations defined by a \code{SpatVector} object. The extraction is
#' performed server-side in GEE and the results are returned as a
#' \code{data.table}. For large point sets, the extraction is automatically
#' chunked to avoid exceeding GEE server memory limits.
#'
#' \strong{Authentication:} This function requires an active Google Earth
#' Engine session. Authenticate using \code{ICESat2VegR::ee_initialize()}
#' before calling this function.
#'
#' @param stack A Google Earth Engine \code{ee$Image} object, typically
#'   the output of \code{\link{ee_build_AlphaEarth_embedding_terrain_stack}}
#'   or \code{\link{ee_build_hls_s1c_terrain_stack}}.
#' @param geom A \code{\link[terra]{SpatVector}} object containing the point
#'   locations at which to extract values. Must be in \code{EPSG:4326}
#'   (WGS84 geographic coordinates). Typically created with
#'   \code{\link[terra]{vect}}.
#' @param scale numeric. The spatial resolution in meters at which to perform
#'   the extraction. Should match the native resolution of the \code{stack}
#'   (e.g., \code{30} for HLS or Landsat-based stacks, \code{10} for
#'   Sentinel-2-based stacks). Default is \code{10}.
#' @param chunk_size numeric. Number of points to process per GEE request.
#'   When the number of points exceeds this value, the extraction is split
#'   into chunks to avoid exceeding GEE server memory limits.
#'   Default is \code{1000}.
#'
#' @return A \code{\link[data.table]{data.table}} with one row per input
#'   point and one column per band in \code{stack}, containing the extracted
#'   pixel values at each point location.
#'
#' @details
#' The function converts the \code{SpatVector} points to GEE geometries
#' and uses \code{ee$Image$sampleRegions()} to extract pixel values at
#' each point. When the number of points exceeds \code{chunk_size}, the
#' points are split into batches and processed sequentially to avoid
#' memory errors on the GEE server.
#'
#' \strong{Authentication:} requires an active GEE session via
#' \code{ICESat2VegR::ee_initialize()}.
#'
#' @examples
#' \dontrun{
#'   Sys.setenv(EE_PROJECT = "your-ee-project-id")
#'   ICESat2VegR::ee_initialize()
#'
#'   library(sf)
#'   library(terra)
#'
#'   # ── AOI from extdata ────────────────────────────────────────
#'   aoi_path <- system.file("extdata",
#'     "aoi_4326.geojson",
#'     package = "ICESat2VegR"
#'   )
#'   boundary <- sf::st_read(aoi_path, quiet = TRUE)
#'   aoi      <- sf::st_as_sfc(sf::st_bbox(boundary))
#'   ee_geom  <- ICESat2VegR:::.as_ee_geom(aoi)
#'
#'   # ── Build embedding + terrain stack ─────────────────────────
#'   stack <- ee_build_AlphaEarth_embedding_terrain_stack(
#'     geom       = ee_geom,
#'     start_year = 2025,
#'     end_year   = 2025
#'   )
#'
#'   # ── Load and sample ICESat-2 segment points ──────────────────
#'   seg_path <- system.file("extdata",
#'     "ATL03_ATL08_example_segments.geojson",
#'     package = "ICESat2VegR"
#'   )
#'   data_raw <- sf::read_sf(seg_path, quiet = TRUE)
#'
#'   df_dt <- data.table::as.data.table(sf::st_drop_geometry(data_raw))
#'   class(df_dt) <- c("icesat2.atl03_atl08_seg_dt", "data.table", "data.frame")
#'
#'   set.seed(1)
#'   df_sampled <- ICESat2VegR::sample(df_dt, method = randomSampling(500))
#'
#'   df_vect <- terra::vect(
#'     as.data.frame(df_sampled),
#'     geom = c("longitude", "latitude"),
#'     crs  = "EPSG:4326"
#'   )
#'
#'   # ── Extract ancillary values at segment locations ────────────
#'   result <- seg_ancillary_extract(
#'     stack      = stack,
#'     geom       = df_vect,
#'     scale      = 30,
#'     chunk_size = 1000
#'   )
#'
#'   head(result)
#'   nrow(result)
#'   names(result)
#' }
#'
#' @import data.table
#' @export
seg_ancillary_extract <- function(stack, geom, scale = 10, chunk_size = 1000) {
  n <- nrow(geom)
  dts <- list()
  k <- 0L

  for (ii in seq.int(1L, n, by = chunk_size)) {
    tail <- min(ii + chunk_size - 1L, n)
    message(sprintf("Processing %d-%d of %d", ii, tail, n))

    sampled <- extract(stack, geom[ii:tail], scale)

    sz <- tryCatch(sampled$size()$getInfo(), error = function(e) 0L)
    if (is.null(sz) || sz == 0L) next

    k <- k + 1L
    dts[[k]] <- ee_to_dt(sampled)
  }

  if (length(dts) == 0L) return(data.table::data.table())
  data.table::rbindlist(dts, use.names = TRUE, fill = TRUE)
}


#' Extract stack values at point locations from Google Earth Engine
#'
#' @description
#' Extracts pixel values from a Google Earth Engine image stack at a set of
#' point locations defined by a \code{SpatVector} object. The extraction is
#' performed server-side in GEE using \code{ee$Image$sampleRegions()} and
#' returns an \code{ee.FeatureCollection} that can be converted to a
#' \code{data.table} using \code{getInfo()}.
#'
#' \strong{Authentication:} This function requires an active Google Earth
#' Engine session. Authenticate using \code{ICESat2VegR::ee_initialize()}
#' before calling this function.
#'
#' @param stack A Google Earth Engine \code{ee$Image} object or a list of
#'   \code{ee$Image} objects. Typically the output of
#'   \code{\link{ee_build_AlphaEarth_embedding_terrain_stack}} or
#'   \code{\link{ee_build_hls_s1c_terrain_stack}}.
#' @param geom A \code{\link[terra]{SpatVector}} object containing the point
#'   locations at which to extract values. Must be in \code{EPSG:4326}
#'   (WGS84 geographic coordinates). Typically created with
#'   \code{\link[terra]{vect}}.
#' @param scale numeric. The spatial resolution in meters at which to perform
#'   the extraction. Should match the native resolution of the \code{stack}
#'   (e.g., \code{30} for HLS or Landsat-based stacks, \code{10} for
#'   Sentinel-2-based stacks).
#'
#' @return An \code{ee.FeatureCollection} with the properties extracted from
#'   the \code{stack} at each point location. Convert to a
#'   \code{\link[data.table]{data.table}} using \code{result$getInfo()}.
#'
#' @details
#' The function internally converts the \code{SpatVector} points to a GEE
#' \code{FeatureCollection} via \code{.points_to_ee_fc()} and composes the
#' input images into a single \code{ee$Image} via \code{.compose_ee_image()}.
#' Extraction is performed with \code{tileScale = 16} to handle large images
#' without running out of GEE memory.
#'
#' For large point sets it is recommended to use
#' \code{\link{seg_ancillary_extract}} instead, which automatically chunks
#' the extraction into batches of \code{chunk_size} points.
#'
#' @examples
#' \dontrun{
#'   Sys.setenv(EE_PROJECT = "your-ee-project-id")
#'   ICESat2VegR::ee_initialize()
#'
#'   library(sf)
#'   library(terra)
#'
#'   # ── AOI and embedding stack ──────────────────────────────────
#'   aoi_path <- system.file("extdata",
#'     "aoi_4326.geojson",
#'     package = "ICESat2VegR"
#'   )
#'   boundary <- sf::st_read(aoi_path, quiet = TRUE)
#'   aoi      <- sf::st_as_sfc(sf::st_bbox(boundary))
#'   ee_geom  <- ICESat2VegR:::.as_ee_geom(aoi)
#'
#'   stack <- ee_build_AlphaEarth_embedding_terrain_stack(
#'     geom       = ee_geom,
#'     start_year = 2025,
#'     end_year   = 2025
#'   )
#'
#'   # ── Load and sample ICESat-2 segment points ──────────────────
#'   seg_path <- system.file("extdata",
#'     "ATL03_ATL08_example_segments.geojson",
#'     package = "ICESat2VegR"
#'   )
#'   data_raw <- sf::read_sf(seg_path, quiet = TRUE)
#'
#'   df_dt <- data.table::as.data.table(sf::st_drop_geometry(data_raw))
#'   class(df_dt) <- c("icesat2.atl03_atl08_seg_dt", "data.table", "data.frame")
#'
#'   set.seed(1)
#'   df_sampled <- ICESat2VegR::sample(df_dt, method = randomSampling(500))
#'
#'   df_vect <- terra::vect(
#'     as.data.frame(df_sampled),
#'     geom = c("longitude", "latitude"),
#'     crs  = "EPSG:4326"
#'   )
#'
#'   # ── Extract values at point locations ────────────────────────
#'   result_ee <- extract(
#'     stack = stack,
#'     geom  = df_vect,
#'     scale = 30
#'   )
#'
#'   # result_ee is an ee.FeatureCollection
#'   class(result_ee)
#'
#'   # Convert to data.table
#'   result_dt <- result_ee$getInfo()
#'   head(result_dt)
#' }
#'
#' @export
extract <- function(stack, geom, scale) {
  img <- .compose_ee_image(stack)    # ee.Image
  fc  <- .points_to_ee_fc(geom)      # ee.FeatureCollection of points

  img$sampleRegions(
    collection = fc,
    scale      = as.integer(scale),
    tileScale  = 16
  )
}
# extract <- function(stack, geom, scale) {
#   final <- c(stack)
#   tempjson <- tempfile(fileext = ".geojson")
#   terra::writeVector(geom, tempjson, filetype = "geojson")
#   parsed <- jsonlite::parse_json(readLines(tempjson))
#   geojson <- ee$FeatureCollection(parsed)

#   sampled <- final$sampleRegions(
#     collection = geojson,
#     scale = scale,
#     tileScale = 16
#   )

#   return(sampled)
# }

# Compose a single ee.Image from common inputs
.compose_ee_image <- function(stack) {
  ee <- reticulate::import("ee", delay_load = FALSE)

  # Already a Python EE object?
  if (inherits(stack, "python.builtin.object")) {
    cls <- try(stack$`__class__`$`__name__`, silent = TRUE)
    if (!inherits(cls, "try-error") && grepl("ImageCollection", cls, ignore.case = TRUE)) {
      return(stack$median())  # reduce collection (change reducer if needed)
    }
    return(stack)  # assume ee.Image
  }

  # List of ee.Image objects
  if (is.list(stack)) {
    imgs <- Filter(function(z) inherits(z, "python.builtin.object"), stack)
    if (!length(imgs)) stop("No ee.Image objects supplied in 'stack'.")
    img <- imgs[[1L]]
    if (length(imgs) > 1L) {
      for (j in 2:length(imgs)) img <- img$addBands(imgs[[j]])
    }
    return(img)
  }

  stop("Unsupported 'stack' type. Provide an ee.Image, an ImageCollection, or a list of ee.Image objects.")
}

# Turn sf / SpatVector points into ee.FeatureCollection
.points_to_ee_fc <- function(geom) {
  ee <- reticulate::import("ee", delay_load = FALSE)

  if (inherits(geom, "sf")) {
    tf <- tempfile(fileext = ".geojson")
    suppressMessages(sf::st_write(geom, tf, quiet = TRUE))
    return(ee$FeatureCollection(jsonlite::read_json(tf, simplifyVector = FALSE)))
  }
  if (inherits(geom, "SpatVector")) {
    tf <- tempfile(fileext = ".geojson")
    terra::writeVector(geom, tf, filetype = "geojson")
    return(ee$FeatureCollection(jsonlite::read_json(tf, simplifyVector = FALSE)))
  }
  stop("geom must be sf or terra::SpatVector (POINT geometries) for sampling.")
}
