#' @export
"[[<-.ee.ee_list.List" <- function(x, i, j, value) {
  x$add(value)
}


#' @export
`mean.ee.ee_list.List` <- function(x, ...) {
  if (!is.null(ee)) {
    return(x$reduce(ee$Reducer$mean()))
  }
}

#' @export
`max.ee.ee_list.List` <- function(x, ..., na.rm = FALSE) {
  args <- list(...)
  stopifnot("Expected only one argument for max for ee.List" = length(args) == 0)


  if (!is.null(ee)) {
    return(x$reduce(ee$Reducer$max()))
  }
}


#' @export
`min.ee.ee_list.List` <- function(x, ..., na.rm = FALSE) {
  args <- list(...)
  stopifnot("Expected only one argument for max for ee.List" = length(args) == 0)

  if (!is.null(ee)) {
    return(x$reduce(ee$Reducer$min()))
  }
}



#' Convert R geometry objects to an Earth Engine geometry
#'
#' @description
#' `.as_ee_geom()` converts a variety of R geometry representations
#' (including `sf`, `sfc`, `terra` objects, bounding boxes, WKT/GeoJSON,
#' and even existing Earth Engine python objects) into a Python
#' `ee$Geometry` suitable for use in Earth Engine workflows.
#'
#' This helper is designed to be flexible and permissive: it accepts
#' most common spatial formats used in R and normalizes them into a
#' standard Earth Engine geometry. When a buffer distance is supplied,
#' the geometry is optionally buffered in meters on the Earth Engine side.
#'
#' @param geom Geometry input. Supported types include:
#'   \itemize{
#'     \item An Earth Engine python object (e.g. `ee$Geometry`, `ee$Feature`,
#'           `ee$FeatureCollection`). Features/FeatureCollections are
#'           reduced to their geometry via `$geometry()`.
#'     \item An `sf` or `sfc` object.
#'     \item A `terra::SpatVector`.
#'     \item A `terra::SpatExtent` (or numeric vector of length 4 specifying
#'           a bounding box as `c(xmin, ymin, xmax, ymax)`).
#'     \item A `data.frame` with longitude/latitude columns (`xcol`, `ycol`).
#'     \item A single-character WKT string.
#'     \item A single-character GeoJSON string.
#'     \item A parsed GeoJSON list with a non-`NULL` `type` element.
#'   }
#' @param xcol Character. Name of the longitude column when `geom` is a
#'   `data.frame`. Default is `"lon"`.
#' @param ycol Character. Name of the latitude column when `geom` is a
#'   `data.frame`. Default is `"lat"`.
#' @param crs Coordinate reference system of the input geometry when `geom`
#'   is an `sf`/`sfc` object or a `data.frame`. Can be any `sf`-compatible
#'   CRS specification (default is EPSG 4326).
#' @param buffer_m Optional numeric. Buffer distance in meters applied to
#'   the resulting Earth Engine geometry. If `NULL` (default) no buffering
#'   is applied. For sf/data.frame inputs, buffering is done on the EE side
#'   using `ee$Geometry$buffer()`.
#'
#' @details
#' For `sf`/`sfc` and `data.frame` inputs, geometries are first transformed
#' to EPSG:4326 and written to a temporary GeoJSON file, which is then read
#' and wrapped into an `ee$FeatureCollection(... )$geometry()`. For
#' `terra::SpatVector` and `terra::SpatExtent` inputs, conversion proceeds
#' via `terra::writeVector()` / `terra::as.polygons()` and an EE rectangle
#' geometry, respectively.
#'
#' Existing Earth Engine python objects are returned as-is, except that
#' `ee$Feature` and `ee$FeatureCollection` objects are coerced to their
#' underlying geometry via `$geometry()`.
#'
#' @return
#' A Python `ee$Geometry` object (or compatible geometry-like EE object)
#' suitable for use in Earth Engine operations.
#'
#' @examples
#' \dontrun{
#'   ee <- reticulate::import("ee", delay_load = FALSE)
#'
#'   # 1) From sf polygon
#'   library(sf)
#'   poly <- st_as_sfc(st_bbox(c(
#'     xmin = -82.4, xmax = -82.2,
#'     ymin =  29.6, ymax =  29.8
#'   ), crs = 4326))
#'
#'   ee_geom1 <- .as_ee_geom(poly)
#'
#'
#'   # 2) From terra::SpatVector
#'   library(terra)
#'   v <- vect(system.file("extdata", "all_boundary.shp", package = "ICESat2VegR"))
#'   ee_geom2 <- .as_ee_geom(v, buffer_m = 30)
#'
#'
#'   # 3) From numeric extent (xmin, ymin, xmax, ymax)
#'   bbox_vec <- c(-82.4, 29.6, -82.2, 29.8)
#'   ee_geom3 <- .as_ee_geom(bbox_vec)
#'
#'
#'   # 4) From data.frame of points
#'   df <- data.frame(
#'     lon = c(-82.3, -82.25),
#'     lat = c(29.65, 29.7)
#'   )
#'   ee_geom4 <- .as_ee_geom(df, buffer_m = 1000)
#' }
#'
#' @export
.as_ee_geom <- function(geom, xcol = "lon", ycol = "lat", crs = 4326, buffer_m = NULL) {
  ee <- reticulate::import("ee", delay_load = FALSE)

  # Already a Python EE object?
  if (inherits(geom, "python.builtin.object")) {
    # Normalize Feature/FeatureCollection → Geometry if needed
    try({
      cls <- geom$`__class__`$`__name__`
      if (grepl("FeatureCollection|Feature", cls, ignore.case = TRUE)) {
        geom <- geom$geometry()
      }
    }, silent = TRUE)
    return(geom)
  }

  # sf / sfc
  if (inherits(geom, "sf") || inherits(geom, "sfc")) {
    gsf <- if (inherits(geom, "sf")) geom else sf::st_sf(geometry = geom)
    if (is.na(sf::st_crs(gsf))) sf::st_crs(gsf) <- crs
    gsf <- sf::st_transform(gsf, 4326)
    tf <- tempfile(fileext = ".geojson")
    suppressMessages(sf::st_write(gsf, tf, quiet = TRUE))
    parsed <- jsonlite::read_json(tf, simplifyVector = FALSE)
    ee_g <- ee$FeatureCollection(parsed)$geometry()
    if (!is.null(buffer_m) && is.finite(buffer_m) && buffer_m > 0) ee_g <- ee_g$buffer(buffer_m)
    return(ee_g)
  }

  # terra::SpatVector
  if (inherits(geom, "SpatVector")) {
    tf <- tempfile(fileext = ".geojson")
    terra::writeVector(geom, tf, filetype = "geojson")
    parsed <- jsonlite::read_json(tf, simplifyVector = FALSE)
    return(ee$FeatureCollection(parsed)$geometry())
  }

  # terra::SpatExtent or numeric bbox c(xmin, ymin, xmax, ymax)
  if (inherits(geom, "SpatExtent") || (is.numeric(geom) && length(geom) == 4L)) {
    bb <- if (inherits(geom, "SpatExtent"))
      c(geom$xmin, geom$ymin, geom$xmax, geom$ymax) else as.numeric(geom)
    return(ee$Geometry$Rectangle(list(bb[1], bb[2], bb[3], bb[4])))
  }

  # data.frame with lon/lat
  if (is.data.frame(geom) && all(c(xcol, ycol) %in% names(geom))) {
    pts <- sf::st_as_sf(geom, coords = c(xcol, ycol), crs = crs) |>
      sf::st_transform(4326)
    tf <- tempfile(fileext = ".geojson")
    suppressMessages(sf::st_write(pts, tf, quiet = TRUE))
    parsed <- jsonlite::read_json(tf, simplifyVector = FALSE)
    ee_g <- ee$FeatureCollection(parsed)$geometry()
    if (!is.null(buffer_m) && is.finite(buffer_m) && buffer_m > 0) ee_g <- ee_g$buffer(buffer_m)
    return(ee_g)
  }

  # WKT string
  if (is.character(geom) && length(geom) == 1L && grepl("^[A-Z]+\\s*\\(", geom)) {
    g <- sf::st_as_sfc(geom, crs = crs)
    return(.as_ee_geom(g, crs = crs, buffer_m = buffer_m))
  }

  # GeoJSON string
  if (is.character(geom) && length(geom) == 1L && grepl("\\{\\s*\"type\"\\s*:", geom)) {
    parsed <- jsonlite::fromJSON(geom, simplifyVector = FALSE)
    return(ee$FeatureCollection(parsed)$geometry())
  }

  # Parsed GeoJSON list
  if (is.list(geom) && !is.null(geom$type)) {
    return(ee$FeatureCollection(geom)$geometry())
  }

  stop("Unsupported 'geom' type. Provide sf/sfc/SpatVector/Extent, bbox, lon/lat data.frame, WKT/GeoJSON, or an EE object.")
}


#' Check if Google Earth Engine is initialized
#'
#' @description
#' `.ee_ping()` verifies that the Earth Engine Python API is initialized and
#' responsive. It attempts to evaluate a trivial Earth Engine expression
#' (`ee$Image$constant(1)$getInfo()`) and returns invisibly if successful.
#' If the call fails, an error is raised with a hint to authenticate and
#' initialize Earth Engine.
#'
#' @param ee The Earth Engine Python module as returned by
#'   `reticulate::import("ee")`.
#'
#' @return
#' Invisibly returns `TRUE` if Earth Engine is initialized and responsive.
#' Otherwise, an error is thrown indicating that Earth Engine must be
#' authenticated and initialized.
#'
#' @examples
#' \dontrun{
#'   library(reticulate)
#'   ee <- import("ee", delay_load = FALSE)
#'
#'   # Will error if EE is not authenticated/initialized
#'   .ee_ping(ee)
#' }
#'
#' @export
.ee_ping <- function(ee) {
  ok <- tryCatch({
    ee$Image$constant(1)$getInfo()
    TRUE
  }, error = function(e) FALSE)

  if (!ok) {
    stop(
      "Earth Engine not initialized. ",
      "Call ee$Authenticate(); ee$Initialize(project = 'your-project')."
    )
  }

  invisible(TRUE)
}
