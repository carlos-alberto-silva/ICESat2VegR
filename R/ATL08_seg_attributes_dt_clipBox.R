# =============================================================================
# Clip ATL08 terrain and canopy attributes by bounding extent
# =============================================================================

#' Clip ATL08 terrain and canopy attributes by bounding extent
#'
#' @description
#' Clips ATL08 terrain and canopy segment attributes within a given bounding
#' extent, defined by a single \code{clip_obj} argument. The clipping extent can
#' be provided as:
#'
#' \itemize{
#'   \item (i) a numeric bounding box, or
#'   \item (ii) a spatial extent object.
#' }
#'
#' @param atl08_seg_att_dt An \code{icesat2.atl08_dt} object (output of
#'   [ICESat2VegR::ATL08_seg_attributes_dt()] function).
#'
#' @param clip_obj Bounding extent used to perform the clipping. Supported
#'   inputs:
#'   \itemize{
#'     \item Numeric vector of length 4: \code{c(xmin, ymin, xmax, ymax)} in
#'           decimal degrees.
#'     \item \code{SpatExtent} (package **terra**): the extent is used directly.
#'   }
#'
#' @return
#' Returns an S4 object of class
#'   [`ICESat2VegR::icesat2.atl08_dt-class`]
#' containing the clipped ATL08 terrain and canopy attributes.
#'
#' @details
#' When \code{clip_obj} is a \code{SpatExtent}, the package \pkg{terra} must be
#' installed. If not found, the function stops with an informative message.
#'
#' @examples
#' \dontrun{
#' # Specifying the path to the ATL08 file
#' atl08_path <- system.file("extdata",
#'   "atl08_clip.h5",
#'   package = "ICESat2VegR"
#' )
#'
#' # Reading ATL08 data (h5 file)
#' atl08_h5 <- ATL08_read(atl08_path)
#'
#' # Extracting ATL08-derived segment attributes
#' atl08_seg_att_dt <- ATL08_seg_attributes_dt(atl08_h5 = atl08_h5)
#'
#' # 1) Using a numeric bounding box: xmin ymin xmax ymax
#' bbox <- c(-103.7604, 59.4672, -103.7600, 59.4680)
#'
#' atl08_seg_att_dt_clip <- ATL08_seg_attributes_dt_clipBox(
#'   atl08_seg_att_dt = atl08_seg_att_dt,
#'   clip_obj = bbox
#' )
#'
#' # 2) Using a SpatExtent
#' # library(terra)
#' # ext <- terra::ext(-103.7604, -103.7600, 59.4672, 59.4680)
#' # atl08_seg_att_dt_clip_ext <- ATL08_seg_attributes_dt_clipBox(
#' #   atl08_seg_att_dt = atl08_seg_att_dt,
#' #   clip_obj = ext
#' # )
#'
#' close(atl08_h5)
#' }
#'
#' @import hdf5r stats
#' @export
ATL08_seg_attributes_dt_clipBox <- function(atl08_seg_att_dt,
                                            clip_obj) {

  if (!inherits(atl08_seg_att_dt, "icesat2.atl08_dt")) {
    stop("atl08_seg_att_dt must be an object of class 'icesat2.atl08_dt'.")
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
  if (any(is.na(atl08_seg_att_dt))) {
    atl08_seg_att_dt <- na.omit(atl08_seg_att_dt)
  }

  mask <-
    atl08_seg_att_dt$longitude >= xmin &
    atl08_seg_att_dt$longitude <= xmax &
    atl08_seg_att_dt$latitude  >= ymin &
    atl08_seg_att_dt$latitude  <= ymax

  mask[!stats::complete.cases(mask)] <- FALSE
  mask <- which(mask)

  newFile <- atl08_seg_att_dt[mask, ]

  # Preserve class
  prepend_class(newFile, c("icesat2.atl08_dt","data.table","data.frame"))
  class(newFile)

  return(newFile)
}
