# --- internal: normalize clip_obj to terra::SpatVector (clip_obj/Multiclip_obj) ---
.as_spatvector_poly_old <- function(clip_obj) {
  if (inherits(clip_obj, "SpatVector")) {
    sv <- clip_obj
  } else if (inherits(clip_obj, c("sf", "sfc"))) {
    if (!requireNamespace("sf", quietly = TRUE)) {
      stop("Package 'sf' is required to pass sf/sfc clip_objs.")
    }
    if (!inherits(clip_obj, "sf")) clip_obj <- sf::st_as_sf(clip_obj)
    gtypes <- unique(as.character(sf::st_geometry_type(clip_obj, by_geometry = TRUE)))
    if (!any(grepl("clip_obj", gtypes, fixed = TRUE))) {
      stop("clip_obj must contain clip_obj or MULTIclip_obj geometries (got: ",
           paste(gtypes, collapse = ", "), ").")
    }
    if (any(!sf::st_is_valid(clip_obj))) clip_obj <- sf::st_make_valid(clip_obj)
    sv <- terra::vect(clip_obj)  # convert to SpatVector (keeps attributes)
  } else {
    stop("`clip_obj` must be a terra::SpatVector or an sf/sfc object.")
  }
  sv <- terra::makeValid(sv)
  gt <- terra::geomtype(sv)
  if (!gt %in% c("clip_objs", "multiclip_objs")) {
    stop("`clip_obj` must be clip_objs/multiclip_objs (got: ", gt, ").")
  }
  sv
}
.as_spatvector_poly <- function(clip_obj) {
  if (inherits(clip_obj, "SpatVector")) {
    sv <- clip_obj
  } else if (inherits(clip_obj, c("sf", "sfc"))) {
    if (!inherits(clip_obj, "sf")) clip_obj <- sf::st_as_sf(clip_obj)
    if (any(!sf::st_is_valid(clip_obj))) clip_obj <- sf::st_make_valid(clip_obj)
    sv <- terra::vect(clip_obj)
  } else {
    stop("`clip_obj` must be a terra::SpatVector or an sf/sfc object.")
  }
  terra::makeValid(sv)
}
#' Clip joined ATL03/ATL08 segment attributes by bounding extent
#'
#' @description
#' Clips an \code{icesat2.atl08_dt} object to a given spatial extent. Only
#' rows whose \code{longitude} and \code{latitude} fall inside the specified
#' bounding extent are retained. This is the bounding-box counterpart to
#' geometry-based clipping via
#' \code{\link{ATL03_ATL08_seg_attributes_dt_clipGeometry}}.
#'
#' @param atl03_atl08_dt An object of class
#'   [`ICESat2VegR::icesat2.atl08_dt-class`], typically produced by
#'   [ATL08_seg_attributes_dt()].
#'
#' @param clip_obj Bounding extent used for clipping. Supported inputs:
#'   \itemize{
#'     \item A numeric vector of length 4:
#'           \code{c(xmin, ymin, xmax, ymax)} in decimal degrees.
#'     \item A \code{\link[terra]{SpatExtent}} object.
#'   }
#'
#' @return An object of the same class as \code{atl03_atl08_dt}, containing
#'   only segments whose \code{longitude} and \code{latitude} fall inside
#'   the bounding extent.
#'
#' @details
#' The function:
#' \enumerate{
#'   \item Validates that the input has \code{longitude} and \code{latitude}
#'         columns.
#'   \item Converts \code{clip_obj} to numeric bounds \code{xmin},
#'         \code{ymin}, \code{xmax}, \code{ymax}.
#'   \item Drops rows with missing coordinates.
#'   \item Returns a subset of \code{atl03_atl08_dt} where
#'         \code{longitude} and \code{latitude} fall inside the bounding box.
#' }
#'
#' @examples
#' \dontrun{
#' # Specifying the path to ATL08 file
#' atl08_path <- system.file("extdata",
#'   "atl08_clip.h5",
#'   package = "ICESat2VegR"
#' )
#'
#' # Reading ATL08 data (h5 file)
#' atl08_h5 <- ATL08_read(atl08_path = atl08_path)
#'
#' # Extracting ATL08 segment attributes
#' atl08_seg_dt <- ATL08_seg_attributes_dt(atl08_h5)
#' head(atl08_seg_dt)
#'
#' # Example 1: Clip using a numeric bounding box (xmin, ymin, xmax, ymax)
#' bbox <- c(-106.571, 41.532, -106.570, 41.537)
#' atl08_seg_clip <- ATL03_ATL08_seg_attributes_dt_clipBox(
#'   atl08_seg_dt,
#'   clip_obj = bbox
#' )
#' head(atl08_seg_clip)
#' nrow(atl08_seg_clip)
#'
#' # Example 2: Clip using a terra SpatExtent object
#' library(terra)
#' ext_obj <- terra::ext(bbox)
#' atl08_seg_clip2 <- ATL03_ATL08_seg_attributes_dt_clipBox(
#'   atl08_seg_dt,
#'   clip_obj = ext_obj
#' )
#' head(atl08_seg_clip2)
#' nrow(atl08_seg_clip2)
#'
#' close(atl08_h5)
#'}
#' @import data.table
#' @export
ATL03_ATL08_seg_attributes_dt_clipBox <- function(atl03_atl08_dt,
                                                  clip_obj) {
  #if (!inherits(atl03_atl08_dt, "icesat2.atl03_atl08_seg_dt")) {
  #  stop("`atl03_atl08_dt` must be of class 'icesat2.atl03_atl08_seg_dt'.")
  #}
  # Accept both classes
  if (!inherits(atl03_atl08_dt, c("icesat2.atl08_dt", "icesat2.atl03_atl08_seg_dt"))) {
    stop("`atl03_atl08_dt` must be of class 'icesat2.atl08_dt' or 'icesat2.atl03_atl08_seg_dt'.")
  }
  if (missing(clip_obj) || is.null(clip_obj)) {
    stop("Argument 'clip_obj' must be provided.")
  }

  # ---------------------------------------------------------------------------
  # Determine bounding box from clip_obj
  # ---------------------------------------------------------------------------

  if (is.numeric(clip_obj)) {
    if (length(clip_obj) != 4L) {
      stop("When numeric, 'clip_obj' must be a vector of length 4: c(xmin, ymin, xmax, ymax).")
    }
    xmin <- clip_obj[1]
    ymin <- clip_obj[2]
    xmax <- clip_obj[3]
    ymax <- clip_obj[4]

  } else if (inherits(clip_obj, "SpatExtent")) {
    if (!requireNamespace("terra", quietly = TRUE)) {
      stop("Package 'terra' is required when using a 'SpatExtent' clip_obj.")
    }
    xmin <- terra::xmin(clip_obj)
    ymin <- terra::ymin(clip_obj)
    xmax <- terra::xmax(clip_obj)
    ymax <- terra::ymax(clip_obj)

  } else {
    stop(
      "Unsupported 'clip_obj' type. ",
      "Supported types are: numeric bbox (length 4) or 'SpatExtent'."
    )
  }

  # ---------------------------------------------------------------------------
  # Input sanity checks
  # ---------------------------------------------------------------------------

  if (!all(c("longitude", "latitude") %in% names(atl03_atl08_dt))) {
    stop("Input must contain 'longitude' and 'latitude' columns.")
  }

  # drop rows with NA lon/lat
  if (anyNA(atl03_atl08_dt$longitude) || anyNA(atl03_atl08_dt$latitude)) {
    keep <- stats::complete.cases(atl03_atl08_dt$longitude,
                                  atl03_atl08_dt$latitude)
    atl03_atl08_dt <- atl03_atl08_dt[keep, ]
  }

  if (nrow(atl03_atl08_dt) == 0L) {
    out <- data.table::data.table()
    class(out) <- class(atl03_atl08_dt)
    return(out)
  }

  # ---------------------------------------------------------------------------
  # Clip by bbox
  # ---------------------------------------------------------------------------

  mask <- atl03_atl08_dt$longitude >= xmin &
    atl03_atl08_dt$longitude <= xmax &
    atl03_atl08_dt$latitude  >= ymin &
    atl03_atl08_dt$latitude  <= ymax

  mask[!stats::complete.cases(mask)] <- FALSE
  idx <- which(mask)

  out <- atl03_atl08_dt[idx, ]
  class(out) <- class(atl03_atl08_dt)

  out
}

#' Clip joined ATL03/ATL08 segment attributes by geometry
#'
#' @description
#' Clips an \code{icesat2.atl08_dt} or \code{icesat2.atl03_atl08_seg_dt}
#' object to the features inside a polygon. Accepts \code{terra::SpatVector}
#' or \code{sf}/\code{sfc} polygons. Performs a fast bounding box pre-clip,
#' then a precise intersection in the polygon CRS.
#'
#' @param atl03_atl08_dt An object of class
#'   [`ICESat2VegR::icesat2.atl08_dt-class`] or
#'   [`ICESat2VegR::icesat2.atl03_atl08_seg_dt-class`],
#'   typically produced by [ATL08_seg_attributes_dt()].
#' @param clip_obj A polygon as \code{terra::SpatVector} or \code{sf}/\code{sfc}.
#' @param split_by character. Optional name of a polygon attribute to copy
#'   into the output as column \code{poly_id}. Default is \code{NULL}.
#'
#' @return An object of the same class as \code{atl03_atl08_dt} containing
#'   only segments that fall inside the polygon. If \code{split_by} is
#'   provided, an additional \code{poly_id} column is added identifying
#'   which polygon feature each segment belongs to.
#'
#' @details
#' The function performs clipping in two steps for efficiency:
#' \enumerate{
#'   \item A fast bounding box pre-clip using
#'     \code{\link{ATL03_ATL08_seg_attributes_dt_clipBox}} to reduce
#'     the number of points before the precise intersection.
#'   \item A precise point-in-polygon intersection using
#'     \code{terra::intersect()}, reprojecting the points to the
#'     polygon CRS if needed.
#' }
#'
#' @examples
#' # Specifying the path to ATL08 file
#' atl08_path <- system.file("extdata",
#'   "atl08_clip.h5",
#'   package = "ICESat2VegR"
#' )
#'
#' # Reading ATL08 data (h5 file)
#' atl08_h5 <- ATL08_read(atl08_path = atl08_path)
#'
#' # Extracting ATL08 segment attributes
#' atl08_seg_dt <- ATL08_seg_attributes_dt(atl08_h5)
#' head(atl08_seg_dt)
#'
#' # Reading polygon from shapefile
#' polygon_filepath <- system.file("extdata",
#'   "clip_geom.shp",
#'   package = "ICESat2VegR"
#' )
#' polygon <- terra::vect(polygon_filepath)
#'
#' # Example 1: Clip by polygon geometry
#' atl08_seg_clip_geom <- ATL03_ATL08_seg_attributes_dt_clipGeometry(
#'   atl08_seg_dt,
#'   clip_obj = polygon
#' )
#' head(atl08_seg_clip_geom)
#' nrow(atl08_seg_clip_geom)
#'
#' # Example 2: Clip by polygon geometry and add polygon id column
#' atl08_seg_clip_geom2 <- ATL03_ATL08_seg_attributes_dt_clipGeometry(
#'   atl08_seg_dt,
#'   clip_obj = polygon,
#'   split_by = "id"
#' )
#' head(atl08_seg_clip_geom2)
#' nrow(atl08_seg_clip_geom2)
#'
#' close(atl08_h5)
#'
#' @import data.table
#' @export
ATL03_ATL08_seg_attributes_dt_clipGeometry <- function(atl03_atl08_dt,
                                                       clip_obj,
                                                       split_by = NULL) {
  if (!inherits(atl03_atl08_dt, c("icesat2.atl08_dt", "icesat2.atl03_atl08_seg_dt"))) {
    stop("`atl03_atl08_dt` must be of class 'icesat2.atl08_dt' or 'icesat2.atl03_atl08_seg_dt'.")
  }
  # check required columns
  if (!all(c("longitude", "latitude") %in% names(atl03_atl08_dt))) {
    stop("Input must contain 'longitude' and 'latitude' columns.")
  }
  # accept terra or sf/sfc, ensure clip_objal, valid
  clip_obj <- .as_spatvector_poly(clip_obj)
  if (!is.null(split_by) && !split_by %in% names(clip_obj)) {
    stop("Field '", split_by, "' not found in clip_obj attributes. Available: ",
         paste(names(clip_obj), collapse = ", "))
  }
  # 1) fast bbox pre-clip
  ex <- terra::ext(clip_obj)
  pre <- ATL03_ATL08_seg_attributes_dt_clipBox(
    atl03_atl08_dt,
    clip_obj = ex
  )
  if (nrow(pre) == 0L) {
    warning("clip_obj extent does not overlap the segment data.")
    out <- data.table::data.table()
    class(out) <- class(atl03_atl08_dt)
    return(out)
  }
  # 2) precise intersection
  pre$nid <- seq_len(nrow(pre))
  # Build points in WGS84 and project to clip_obj CRS if needed
  pts <- terra::vect(
    as.data.frame(pre),
    geom = c("longitude", "latitude"),
    crs  = "EPSG:4326"
  )
  if (!terra::same.crs(pts, clip_obj)) {
    pts <- terra::project(pts, terra::crs(clip_obj))
  }
  inter <- terra::intersect(terra::makeValid(pts), clip_obj)
  if (nrow(inter) == 0L) {
    warning("No segments fall inside the clip_obj (after precise intersection).")
    out <- data.table::data.table()
    class(out) <- class(atl03_atl08_dt)
    return(out)
  }
  out <- pre[inter$nid, ]
  if (!is.null(split_by)) out$poly_id <- inter[[split_by]]
  # cleanup & preserve class
  out$nid <- NULL
  class(out) <- class(atl03_atl08_dt)
  out
}
