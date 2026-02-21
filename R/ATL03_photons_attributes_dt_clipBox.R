#' Clip ATL03 photons by bounding extent
#'
#' @description
#' Clips ATL03 photon attributes within a given bounding extent, defined by a
#' single `clip_obj` argument. The clipping extent can be provided as:
#'
#' \itemize{
#'   \item (i) a numeric bounding box, or
#'   \item (ii) a spatial extent object.
#' }
#'
#' @param atl03_photons_dt An `icesat2.atl03_dt` object (output of
#'   [ATL03_photons_attributes_dt()]).
#'
#' @param clip_obj Bounding extent used to perform the clipping. Supported
#'   inputs:
#'   \itemize{
#'     \item Numeric vector of length 4: `c(xmin, ymin, xmax, ymax)` in
#'           decimal degrees.
#'     \item `SpatExtent` (package **terra**): the extent is used directly.
#'   }
#'
#' @return
#' Returns a subset of the original `icesat2.atl03_dt` object containing
#' only photons within the bounding box.
#'
#' @details
#' When `clip_obj` is a `SpatExtent`, the package \pkg{terra} must be
#' installed. If not found, the function stops with an informative message.
#'
#' @seealso
#'  \url{https://icesat-2.gsfc.nasa.gov/sites/default/files/page_files/ICESat2_ATL03_ATBD_r006.pdf}
#'
#' @examples
#' \dontrun{
#' atl03_path <- system.file("extdata", "atl03_clip.h5",
#'                           package = "ICESat2VegR")
#'
#' atl03_h5 <- ATL03_read(atl03_path = atl03_path)
#' atl03_photons_dt <- ATL03_photons_attributes_dt(atl03_h5 = atl03_h5)
#'
#' # 1) Using a numeric bbox: xmin ymin xmax ymax
#' bbox <- c(-106.57, 41.53, -106.5698, 41.54)
#' atl03_clip_bbox <- ATL03_photons_attributes_dt_clipBox(
#'   atl03_photons_dt = atl03_photons_dt,
#'   clip_obj = bbox
#' )
#'
#' # 2) Using a SpatExtent
#' # library(terra)
#' # ext <- terra::ext(-106.57, -106.5698, 41.53, 41.54)
#' # atl03_clip_ext <- ATL03_photons_attributes_dt_clipBox(
#' #   atl03_photons_dt = atl03_photons_dt,
#' #   clip_obj = ext
#' # )
#'
#' close(atl03_h5)
#' }
#'
#' @import hdf5r stats
#' @export
ATL03_photons_attributes_dt_clipBox <- function(atl03_photons_dt,
                                                clip_obj) {

  if (!inherits(atl03_photons_dt, "icesat2.atl03_dt")) {
    stop("atl03_photons_dt must be an object of class 'icesat2.atl03_dt'.")
  }

  if (missing(clip_obj) || is.null(clip_obj)) {
    stop("Argument 'clip_obj' must be provided.")
  }

  # ---------------------------------------------------------------------------
  # Determine bounding box from clip_obj
  # ---------------------------------------------------------------------------

  # Case 1: numeric bbox
  if (is.numeric(clip_obj)) {
    if (length(clip_obj) != 4L) {
      stop("When numeric, 'clip_obj' must be a vector of length 4: c(xmin, ymin, xmax, ymax).")
    }
    xmin <- clip_obj[1]
    ymin <- clip_obj[2]
    xmax <- clip_obj[3]
    ymax <- clip_obj[4]

    # Case 2: terra SpatExtent
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
  # Clean and clip
  # ---------------------------------------------------------------------------
  if (any(is.na(atl03_photons_dt))) {
    atl03_photons_dt <- na.omit(atl03_photons_dt)
  }

  mask <-
    atl03_photons_dt$lon_ph >= xmin &
    atl03_photons_dt$lon_ph <= xmax &
    atl03_photons_dt$lat_ph >= ymin &
    atl03_photons_dt$lat_ph <= ymax

  mask[!stats::complete.cases(mask)] <- FALSE
  mask <- which(mask)

  clipped <- atl03_photons_dt[mask, ]

  return(clipped)
}
