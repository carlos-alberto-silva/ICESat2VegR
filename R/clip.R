# =============================================================================
# Generic: clip
# =============================================================================

#' Clip ICESat-2 data (h5, attributes, or ATL03–ATL08 join) by box or geometry
#'
#' @description
#' Unified clipping interface for ICESat-2 ATL03/ATL08 data. The generic
#' `clip()` dispatches on the class of `x` and internally calls appropriate
#' `_clipBox()` or `_clipGeometry()` helpers.
#'
#' @param x ICESat-2 object (ATL03/ATL08 h5, attributes, or joined table)
#'
#' @param clip_obj Clipping object. Supported for box-based clipping:
#'   * numeric bounding box: c(xmin, ymin, xmax, ymax)
#'   * `terra::SpatExtent`
#'
#'   For geometry-based clipping:
#'   * `sf`, `sfc`, `terra::SpatVector`, etc.
#'
#'   To clip by a spatial object's extent, extract its bbox first.
#'
#' @export
setGeneric(
  "clip",
  function(x, clip_obj, ...) standardGeneric("clip")
)

# =============================================================================
# Methods: clip() for ATL03/ATL08 HDF5 handles
# =============================================================================

# --- ATL03 H5 by numeric bbox ---
setMethod(
  "clip",
  signature = c(x = "icesat2.atl03_h5", clip_obj = "numeric"),
  function(x, clip_obj, ...) {
    ATL03_h5_clipBox(x, clip_obj = clip_obj, ...)
  }
)

# --- ATL03 H5 by SpatExtent ---
setMethod(
  "clip",
  signature = c(x = "icesat2.atl03_h5", clip_obj = "SpatExtent"),
  function(x, clip_obj, ...) {
    ATL03_h5_clipBox(x, clip_obj = clip_obj, ...)
  }
)

# --- ATL03 H5 by geometry ---
setMethod(
  "clip",
  signature = c(x = "icesat2.atl03_h5", clip_obj = "ANY"),
  function(x, clip_obj, ...) {
    ATL03_h5_clipGeometry(x, clip_obj = clip_obj, ...)
  }
)

# --- ATL08 H5 by numeric bbox ---
setMethod(
  "clip",
  signature = c(x = "icesat2.atl08_h5", clip_obj = "numeric"),
  function(x, clip_obj, ...) {
    ATL08_h5_clipBox(x, clip_obj = clip_obj, ...)
  }
)

# --- ATL08 H5 by SpatExtent ---
setMethod(
  "clip",
  signature = c(x = "icesat2.atl08_h5", clip_obj = "SpatExtent"),
  function(x, clip_obj, ...) {
    ATL08_h5_clipBox(x, clip_obj = clip_obj, ...)
  }
)

# --- ATL08 H5 by geometry ---
setMethod(
  "clip",
  signature = c(x = "icesat2.atl08_h5", clip_obj = "ANY"),
  function(x, clip_obj, ...) {
    ATL08_h5_clipGeometry(x, clip_obj = clip_obj, ...)
  }
)

# =============================================================================
# Methods: clip() for ATL03/ATL08 attribute tables
# =============================================================================

# --- ATL03 DT by numeric bbox ---
setMethod(
  "clip",
  signature = c(x = "icesat2.atl03_dt", clip_obj = "numeric"),
  function(x, clip_obj, ...) {
    ATL03_photons_attributes_dt_clipBox(x, clip_obj = clip_obj, ...)
  }
)

# --- ATL03 DT by SpatExtent ---
setMethod(
  "clip",
  signature = c(x = "icesat2.atl03_dt", clip_obj = "SpatExtent"),
  function(x, clip_obj, ...) {
    ATL03_photons_attributes_dt_clipBox(x, clip_obj = clip_obj, ...)
  }
)

# --- ATL03 DT by geometry ---
setMethod(
  "clip",
  signature = c(x = "icesat2.atl03_dt", clip_obj = "ANY"),
  function(x, clip_obj, ...) {
    ATL03_photons_attributes_dt_clipGeometry(x, clip_obj = clip_obj, ...)
  }
)

# --- ATL08 DT by numeric bbox ---
setMethod(
  "clip",
  signature = c(x = "icesat2.atl08_dt", clip_obj = "numeric"),
  function(x, clip_obj, ...) {
    ATL08_seg_attributes_dt_clipBox(x, clip_obj = clip_obj, ...)
  }
)

# --- ATL08 DT by SpatExtent ---
setMethod(
  "clip",
  signature = c(x = "icesat2.atl08_dt", clip_obj = "SpatExtent"),
  function(x, clip_obj, ...) {
    ATL08_seg_attributes_dt_clipBox(x, clip_obj = clip_obj, ...)
  }
)

# --- ATL08 DT by geometry ---
setMethod(
  "clip",
  signature = c(x = "icesat2.atl08_dt", clip_obj = "ANY"),
  function(x, clip_obj, ...) {
    ATL08_seg_attributes_dt_clipGeometry(x, clip_obj = clip_obj, ...)
  }
)

# =============================================================================
# Methods: clip() for joined ATL03–ATL08 photons
# =============================================================================

# --- ATL03-ATL08 joined DT by numeric bbox ---
setMethod(
  "clip",
  signature = c(x = "icesat2.atl03atl08_dt", clip_obj = "numeric"),
  function(x, clip_obj, ...) {
    ATL03_ATL08_photons_attributes_dt_clipBox(
      atl03_atl08_dt = x,
      clip_obj       = clip_obj,
      ...
    )
  }
)

# --- ATL03-ATL08 joined DT by SpatExtent ---
setMethod(
  "clip",
  signature = c(x = "icesat2.atl03atl08_dt", clip_obj = "SpatExtent"),
  function(x, clip_obj, ...) {
    ATL03_ATL08_photons_attributes_dt_clipBox(
      atl03_atl08_dt = x,
      clip_obj       = clip_obj,
      ...
    )
  }
)

# --- ATL03-ATL08 joined DT by geometry ---
setMethod(
  "clip",
  signature = c(x = "icesat2.atl03atl08_dt", clip_obj = "ANY"),
  function(x, clip_obj, ...) {
    ATL03_ATL08_photons_attributes_dt_clipGeometry(
      atl03_atl08_dt = x,
      clip_obj       = clip_obj,
      ...
    )
  }
)
