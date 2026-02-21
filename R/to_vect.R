#' Convert ICESat-2 classes, data.frame/data.table, and sf to terra::SpatVector
#'
#' Generic function to convert ICESat-2 objects, plain tables (data.frame / data.table),
#' and `sf` objects to a [`terra::SpatVector-class`].
#'
#' For plain tables, you can either:
#' 1) Provide `lon=` and `lat=` column names, or
#' 2) Rely on defaults (`longitude` / `latitude`), or
#' 3) Use common ICESat-2 column names (`lon_ph`/`lat_ph`, `reference_photon_lon`/`reference_photon_lat`).
#'
#' @param x The ICESat-2, `data.frame`/`data.table`, or `sf` object to convert.
#' @param ... Additional parameters passed to methods.
#'
#' @return A [`terra::SpatVector-class`] object.
#'
#' @examples
#' # ICESat2VegR examples (package objects keep their classes):
#' atl03_path <- system.file("extdata", "atl03_clip.h5", package = "ICESat2VegR")
#' atl03_h5 <- ATL03_read(atl03_path = atl03_path)
#' atl03_segment_dt <- ATL03_seg_metadata_dt(atl03_h5 = atl03_h5)
#' atl03_segment_vect <- to_vect(atl03_segment_dt)
#' terra::plot(atl03_segment_vect, col = atl03_segment_vect$segment_ph_cnt)
#' close(atl03_h5)
#'
#' # Plain data.frame / data.table usage:
#' # df <- data.frame(longitude = c(-84, -84.1), latitude = c(29.6, 29.7), z = 1:2)
#' # v  <- to_vect(df)
#' # v2 <- to_vect(df, lon = "longitude", lat = "latitude", crs = "EPSG:4326")
#'
#' @include class.icesat2.R
#' @export
setGeneric(
  "to_vect",
  def = function(x, ...) {
    standardGeneric("to_vect")
  }
)

# ---- S4 registration for S3 classes (needed for S4 dispatch) ----
# NOTE: data.table is S3; to enable setMethod("to_vect", "data.table", ...)
# we must register it with setOldClass.
# data.table's class vector is c("data.table","data.frame"), so register both forms.
if (!methods::isClass("data.frame")) {
  methods::setOldClass("data.frame")
}
if (!methods::isClass("data.table")) {
  methods::setOldClass(c("data.table", "data.frame"))
}

# ---- internal helper ---- 
.generic_to_vect <- function(x, lon, lat, crs = "EPSG:4326") {
  if (!requireNamespace("terra", quietly = TRUE)) {
    stop("Package 'terra' is required but not installed.")
  }

  # Convert to data.frame if needed
  if (inherits(x, "data.table")) x <- as.data.frame(x)
  if (!inherits(x, "data.frame")) {
    stop("Expected a data.frame/data.table-like object; got: ", paste(class(x), collapse = ", "))
  }

  # Check coordinate columns
  if (!all(c(lon, lat) %in% names(x))) {
    stop("Columns '", lon, "' and/or '", lat, "' not found in data.")
  }

  suppressWarnings({
    x[[lon]] <- as.numeric(x[[lon]])
    x[[lat]] <- as.numeric(x[[lat]])
  })

  ok <- is.finite(x[[lon]]) & is.finite(x[[lat]])
  x <- x[ok, , drop = FALSE]

  if (nrow(x) == 0) {
    stop("No finite coordinates to export.")
  }

  if (grepl("4326|WGS", toupper(crs))) {
    if (any(
      x[[lon]] < -180 | x[[lon]] > 180 |
      x[[lat]] < -90  | x[[lat]] > 90,
      na.rm = TRUE
    )) {
      warning("Some coordinates are outside WGS84 bounds - check CRS or input values.")
    }
  }

  terra::vect(x, geom = c(lon, lat), crs = crs)
}

# Choose lon/lat columns for plain tables
# Priority:
# 1) user-supplied lon/lat
# 2) common lon/lat
# 3) ICESat-2 typical columns
# 4) fail with a helpful message
.resolve_lonlat <- function(x, lon = NULL, lat = NULL) {
  nms <- names(x)

  if (!is.null(lon) && !is.null(lat)) return(list(lon = lon, lat = lat))

  candidates <- list(
    c("longitude", "latitude"),
    c("lon", "lat"),
    c("x", "y"),
    c("lon_ph", "lat_ph"),
    c("reference_photon_lon", "reference_photon_lat")
  )

  for (pair in candidates) {
    if (all(pair %in% nms)) return(list(lon = pair[1], lat = pair[2]))
  }

  stop(
    "Could not infer coordinate columns. Provide `lon=` and `lat=`.\n",
    "Available columns include: ", paste(utils::head(nms, 30), collapse = ", "),
    if (length(nms) > 30) " ..."
  )
}

# ---- Methods for plain tables ----

#' @rdname to_vect
#' @param lon,lat For `data.frame`/`data.table`, the column names to use as coordinates.
#' @param crs CRS string for the resulting SpatVector. Default is "EPSG:4326".
#' @export
setMethod(
  "to_vect",
  "data.frame",
  function(x, lon = NULL, lat = NULL, crs = "EPSG:4326", ...) {
    ll <- .resolve_lonlat(x, lon = lon, lat = lat)
    .generic_to_vect(x, lon = ll$lon, lat = ll$lat, crs = crs)
  }
)

#' @rdname to_vect
#' @export
setMethod(
  "to_vect",
  "data.table",
  function(x, lon = NULL, lat = NULL, crs = "EPSG:4326", ...) {
    ll <- .resolve_lonlat(x, lon = lon, lat = lat)
    .generic_to_vect(x, lon = ll$lon, lat = ll$lat, crs = crs)
  }
)

# ---- ICESat-2 class methods ----

#' @rdname to_vect
#' @export
setMethod(
  "to_vect",
  "icesat2.atl03_seg_dt",
  function(x, ...) {
    .generic_to_vect(x, "reference_photon_lon", "reference_photon_lat")
  }
)

#' @rdname to_vect
#' @export
setMethod(
  "to_vect",
  "icesat2.atl08_dt",
  function(x, ...) {
    .generic_to_vect(x, "longitude", "latitude")
  }
)

#' @rdname to_vect
#' @export
setMethod(
  "to_vect",
  "icesat2.atl03atl08_dt",
  function(x, ...) {
    .generic_to_vect(x, "lon_ph", "lat_ph")
  }
)

#' @rdname to_vect
#' @export
setMethod(
  "to_vect",
  "icesat2.atl03_dt",
  function(x, ...) {
    .generic_to_vect(x, "lon_ph", "lat_ph")
  }
)

#' @rdname to_vect
#' @export
setMethod(
  "to_vect",
  "icesat2.atl03_atl08_seg_dt",
  function(x, ...) {
    .generic_to_vect(x, "longitude", "latitude")
  }
)

# ---- sf method ----

#' Convert sf to terra::SpatVector
#'
#' @param x An `sf` object (POINT geometry recommended).
#' @param target_crs Optional CRS string (e.g., "EPSG:4326") to reproject the output.
#' @param ... Ignored.
#'
#' @return A `terra::SpatVector`
#' @rdname to_vect
#' @export
setMethod(
  "to_vect",
  "sf",
  function(x, target_crs = NULL, ...) {
    if (!requireNamespace("terra", quietly = TRUE)) {
      stop("Package 'terra' is required but not installed.")
    }

    v <- terra::vect(x)

    if (!is.null(target_crs)) {
      v <- terra::project(v, target_crs)
    }

    v
  }
)
