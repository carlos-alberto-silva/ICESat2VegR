#' Convert vector data to Google Earth Engine FeatureCollection (no rgee)
#'
#' @description
#' `vect_as_ee()` converts vector data to a Google Earth Engine
#' `ee$FeatureCollection` using **reticulate** to call the official Python
#' Earth Engine API-**without** relying on **rgee**. It accepts `sf`/`sfc`/`sfg`
#' objects or [`terra::SpatVector-class`], validates and reprojects geometries, and
#' either returns an in-memory FeatureCollection or exports to an EE Asset.
#'
#' @param x An input vector object: an `sf`, `sfc`, or `sfg` object, or a
#'   [`terra::SpatVector-class`]. For `sfg`/`sfc`, a simple point/line/polygon geometry
#'   is accepted; for `SpatVector`, it will be converted to `sf`.
#' @param via Character; one of `c("getInfo","getInfo_to_asset")`.
#'   * `getInfo` returns an in-memory `ee$FeatureCollection` built from the
#'   GeoJSON of `x`.
#'   * `getInfo_to_asset` creates an EE export task to save the collection to
#'   an Earth Engine **Asset** (requires `assetId`).
#' @param assetId Character; EE asset id (e.g., `users/you/my_fc`). **Required**
#'   when `via = "getInfo_to_asset"`.
#' @param overwrite Logical; if `TRUE`, attempt to delete an existing EE asset at
#'   `assetId` before export.
#' @param proj Character CRS string (e.g., `EPSG:4326`). Input will be
#'   transformed to this CRS before conversion. Default is `EPSG:4326`.
#' @param make_valid Logical; if `TRUE`, attempt to repair invalid geometries
#'   (`sf::st_make_valid()` / `terra::makeValid()`).
#' @param quiet Logical; if `FALSE`, print progress messages (e.g., task state).
#' @param monitoring Logical; when exporting to asset, if `TRUE` poll the task
#'   until completion or timeout, otherwise return immediately after starting.
#' @param poll_interval_sec Numeric; seconds to wait between task status checks
#'   when `monitoring = TRUE`. Default `5`.
#' @param poll_timeout_sec Numeric; maximum seconds to keep polling before
#'   timing out. Default `3600` (1 hour).
#'
#' @return
#' If `via = "getInfo"`, returns an in-memory `ee$FeatureCollection` object
#' (reticulate Python object). If `via = "getInfo_to_asset"`, returns an
#' `ee$FeatureCollection` **referencing** `assetId` (after successful export),
#' or returns immediately if `monitoring = FALSE`.
#'
#' @details
#' * The function converts `x` to `sf`, enforces a defined CRS, optionally fixes
#'   invalid geometries, and transforms to `proj` (default WGS84).
#' * Properties are carried alongside geometry; however, Earth Engine does not
#'   accept `POSIX*` timestamp columns or property names containing `'.'`. The
#'   function stops with a clear error if such columns are found-convert them to
#'   character or rename before calling.
#' * GeoJSON is built in-memory (via **geojsonsf** when available) or via a
#'   temporary `.geojson` file as a fallback, then parsed to a list for
#'   `ee$FeatureCollection`.
#' * Requires a properly initialized EE Python environment (`ee.Initialize()`
#'   in the active Python session used by **reticulate**).
#'
#' @section Errors & Constraints:
#' * Undefined CRS: the function stops if `x` has no CRS set.
#' * Unsupported columns: the function stops if any property is `POSIX*` or if
#'   names contain dots (`.`).
#' * Export failures: if `via = "getInfo_to_asset"` and the EE task fails or is
#'   cancelled, an error is thrown when `monitoring = TRUE`.
#'
#' @examples
#' \dontrun{
#' # -- Prerequisites (Python side):
#' # import ee; ee.Initialize()
#' #
#' # Example with sf POINTS
#' library(sf)
#' pts <- st_as_sf(data.frame(
#'   id = 1:3,
#'   x  = c(-84.02171, -84.02025, -84.02026),
#'   y  = c(31.29891, 31.31945, 31.31938)
#' ), coords = c("x","y"), crs = "EPSG:4326")
#'
#' # In-memory FeatureCollection:
#' fc <- vect_as_ee(pts, via = "getInfo")
#'
#' # Export to an EE asset (monitor until done):
#' # fc_asset <- vect_as_ee(
#' #   pts,
#' #   via = "getInfo_to_asset",
#' #   assetId = "users/you/demo_points",
#' #   overwrite = TRUE,
#' #   monitoring = TRUE
#' # )
#' }
#'
#' @seealso
#' * Earth Engine Python API docs: `ee$FeatureCollection`
#' * Geometry repair: [sf::st_make_valid()], [terra::makeValid()]
#' * GeoJSON helpers: **geojsonsf**, [sf::st_write()]
#'
#' @importFrom sf st_as_sf st_sfc st_sf st_crs st_make_valid st_transform st_write
#' @importFrom terra vect makeValid project
#' @importFrom jsonlite fromJSON
#' @importFrom utils read.table
#' @export
vect_as_ee <- function(
    x,
    via = c("getInfo", "getInfo_to_asset"),
    assetId = NULL,
    overwrite = TRUE,
    proj = "EPSG:4326",
    make_valid = TRUE,
    quiet = FALSE,
    monitoring = TRUE,
    poll_interval_sec = 5,
    poll_timeout_sec = 3600
) {
  via <- match.arg(via)

  # deps
  if (!requireNamespace("sf", quietly = TRUE)) stop("Package 'sf' is required.")

  # Accept terra::SpatVector
  if (inherits(x, "SpatVector")) {
    if (!requireNamespace("terra", quietly = TRUE))
      stop("Input is terra::SpatVector but package 'terra' is not installed.")
    if (isTRUE(make_valid)) x <- tryCatch(terra::makeValid(x), error = function(e) x)
    x <- sf::st_as_sf(x)
  }

  # Normalize to sf without any data.frame coercion
  if (inherits(x, "sfg")) x <- sf::st_sfc(x)
  if (inherits(x, "sfc")) x <- sf::st_sf(geometry = x)
  if (!inherits(x, "sf")) stop("x must be of class sf/sfc/sfg or terra::SpatVector.")

  # CRS & validity
  if (is.na(sf::st_crs(x)$wkt)) stop("Input CRS is undefined; set with sf::st_set_crs().")
  if (isTRUE(make_valid)) x <- tryCatch(sf::st_make_valid(x), error = function(e) x)
  x <- sf::st_transform(x, proj)

  # Forbid POSIX* columns and dots in names (EE constraints)
  props_df <- x
  sf::st_geometry(props_df) <- NULL
  if (ncol(props_df) > 0) {
    bad_time <- vapply(props_df, function(col) any(inherits(col, c("POSIXlt","POSIXct","POSIXt"))), logical(1))
    if (any(bad_time)) {
      stop("POSIX* columns are not supported; convert to character: ",
           paste(names(props_df)[bad_time], collapse = ", "))
    }
    bad_names <- grepl("\\.", names(props_df))
    if (any(bad_names)) {
      stop("Column names with '.' are not supported in Earth Engine: ",
           paste(names(props_df)[bad_names], collapse = ", "))
    }
  }

  # ---- GeoJSON builders (no st_as_geojson dependency) ----
  sf_to_geojson_list <- function(sfobj) {
    # Prefer geojsonsf for in-memory conversion
    if (requireNamespace("geojsonsf", quietly = TRUE)) {
      # Full FeatureCollection as a JSON string
      gj_str <- geojsonsf::sf_geojson(sfobj, simplify = FALSE, digits = 15)
      return(jsonlite::fromJSON(gj_str, simplifyVector = FALSE))
    }

    # Fallback: write to a temporary GeoJSON file, then read it back
    tmp <- tempfile(fileext = ".geojson")
    on.exit(unlink(tmp), add = TRUE)
    sf::st_write(sfobj, tmp, driver = "GeoJSON", quiet = TRUE)
    jsonlite::fromJSON(readLines(tmp, warn = FALSE), simplifyVector = FALSE)
  }

  fc_geojson <- sf_to_geojson_list(x)

  # Earth Engine (Python) via reticulate
  if (!requireNamespace("reticulate", quietly = TRUE)) stop("Package 'reticulate' is required.")
  ee <- reticulate::import("ee", delay_load = FALSE)

  ee_fc_from_geojson <- function(gj) ee$FeatureCollection(gj)

  if (via == "getInfo") {
    # in-memory FeatureCollection
    return(ee_fc_from_geojson(fc_geojson))
  }

  # via == "getInfo_to_asset"
  if (is.null(assetId)) stop("assetId must be provided for 'getInfo_to_asset'.")

  fc <- ee_fc_from_geojson(fc_geojson)

  # overwrite if requested
  if (isTRUE(overwrite)) {
    existing <- try(ee$data$getAsset(assetId), silent = TRUE)
    if (!inherits(existing, "try-error") && !is.null(existing)) {
      if (!quiet) message("Deleting existing asset: ", assetId)
      try(ee$data$deleteAsset(assetId), silent = TRUE)
    }
  }

  desc <- paste0("vect_as_ee_export_", format(Sys.time(), "%Y%m%d_%H%M%S"))
  task <- ee$batch$Export$table$toAsset(
    collection  = fc,
    description = desc,
    assetId     = assetId
  )
  task$start()

  if (!isTRUE(monitoring)) {
    if (!quiet) message("Export started (not monitoring): ", desc)
    return(ee$FeatureCollection(assetId))
  }

  # simple polling
  t0 <- Sys.time()
  `%||%` <- function(a, b) if (is.null(a)) b else a
  repeat {
    st <- try(task$status(), silent = TRUE)
    state <- if (inherits(st, "try-error")) "UNKNOWN" else (st$state %||% "UNKNOWN")
    if (!quiet) message("Task state: ", state)
    if (state %in% c("COMPLETED","FAILED","CANCELLED")) break
    if (difftime(Sys.time(), t0, units = "secs") > poll_timeout_sec) stop("Export monitor timed out.")
    Sys.sleep(poll_interval_sec)
  }
  if (state != "COMPLETED") stop("Export failed: state = ", state)
  ee$FeatureCollection(assetId)
}
