#' Clip joined ATL03 and ATL08 photons by bounding extent
#'
#' @description
#' Clips joined ATL03 and ATL08 photon attributes within a given bounding
#' extent, defined by a single `clip_obj` argument. The clipping extent can
#' be provided as:
#'
#' \itemize{
#'   \item (i) a numeric bounding box, or
#'   \item (ii) a spatial extent object.
#' }
#'
#' @param atl03_atl08_dt An S4 object of class
#'   [`ICESat2VegR::icesat2.atl03atl08_dt-class`] containing joined ATL03 and
#'   ATL08 data (output of
#'   [ICESat2VegR::ATL03_ATL08_photons_attributes_dt_join()] function).
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
#' Returns an S4 object of class
#' [`ICESat2VegR::icesat2.atl03atl08_dt-class`] containing a subset of the
#' joined ATL03 and ATL08 photon attributes restricted to the clipping extent.
#'
#' @seealso \url{https://icesat-2.gsfc.nasa.gov/sites/default/files/page_files/ICESat2_ATL03_ATBD_r006.pdf}
#' @seealso \url{https://icesat-2.gsfc.nasa.gov/sites/default/files/page_files/ICESat2_ATL08_ATBD_r006.pdf}
#'
#' @details
#' When `clip_obj` is a `SpatExtent`, the package \pkg{terra} must be
#' installed. If not found, the function stops with an informative message.
#'
#' @examples
#' \dontrun{
#' # Specifying the path to ATL03 and ATL08 files
#' atl03_path <- system.file("extdata", "atl03_clip.h5",
#'                           package = "ICESat2VegR")
#'
#' atl08_path <- system.file("extdata", "atl08_clip.h5",
#'                           package = "ICESat2VegR")
#'
#' # Reading ATL03 and ATL08 data (h5 files)
#' atl03_h5 <- ATL03_read(atl03_path = atl03_path)
#' atl08_h5 <- ATL08_read(atl08_path = atl08_path)
#'
#' # Joining ATL03 and ATL08 photons and heights
#' atl03_atl08_dt <- ATL03_ATL08_photons_attributes_dt_join(atl03_h5, atl08_h5)
#' head(atl03_atl08_dt)
#'
#' # 1) Using a numeric bounding box: xmin ymin xmax ymax
#' bbox <- c(-103.7604, 59.4672, -103.7600, 59.4680)
#'
#' atl03_atl08_dt_clip <- ATL03_ATL08_photons_attributes_dt_clipBox(
#'   atl03_atl08_dt = atl03_atl08_dt,
#'   clip_obj = bbox
#' )
#' head(atl03_atl08_dt_clip)
#'
#' # 2) Using a SpatExtent (example)
#' # library(terra)
#' # ext <- terra::ext(-103.7604, -103.7600, 59.4672, 59.4680)
#' # atl03_atl08_dt_clip_ext <- ATL03_ATL08_photons_attributes_dt_clipBox(
#' #   atl03_atl08_dt = atl03_atl08_dt,
#' #   clip_obj = ext
#' # )
#'
#' close(atl03_h5)
#' close(atl08_h5)
#' }
#'
#' @import hdf5r stats
#' @export
ATL03_ATL08_photons_attributes_dt_clipBox <- function(atl03_atl08_dt,
                                                      clip_obj) {

  if (!inherits(atl03_atl08_dt, "icesat2.atl03atl08_dt")) {
    stop("atl03_atl08_dt must be an object of class 'icesat2.atl03atl08_dt'.")
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
  if (any(is.na(atl03_atl08_dt))) {
    atl03_atl08_dt <- na.omit(atl03_atl08_dt)
  }

  mask <-
    atl03_atl08_dt$lon_ph >= xmin &
    atl03_atl08_dt$lon_ph <= xmax &
    atl03_atl08_dt$lat_ph >= ymin &
    atl03_atl08_dt$lat_ph <= ymax

  mask[!stats::complete.cases(mask)] <- FALSE
  mask <- which(mask)

  newFile <- atl03_atl08_dt[mask, ]

  # Preserve class
  prepend_class(newFile, "icesat2.atl03atl08_dt")

  return(newFile)
}
