#' Normalize photon heights relative to estimated ground elevation
#'
#' @description
#' Normalizes ATL03 and ATL08 photon heights (\code{ph_h}) by subtracting
#' the estimated ground elevation at each photon location. Ground elevation
#' is estimated internally using \code{\link{ATL03_ATL08_photons_seg_dt_fitground}}.
#' The function updates the \code{ph_h} column in place and returns the
#' modified \code{atl03_atl08_seg_dt} object.
#'
#' @param atl03_atl08_seg_dt An S4 object of class
#'   [`ICESat2VegR::icesat2.atl03_atl08_seg_dt-class`] containing ATL03
#'   and ATL08 data (output of [ATL03_ATL08_segment_create()] function).
#' @param smoothing_window numeric. The smoothing window size in meters
#'   for aggregating ground photons before interpolation.
#'   Default is \code{NA}, which uses the ATBD adaptive window size.
#'   See Details for more information.
#' @param smoothing_func function. The aggregation function applied to
#'   ground photon elevations within each smoothing window.
#'   Default is \code{median}.
#' @param interpolation_func function. The interpolation function used
#'   to estimate ground elevation from the smoothed ground photons.
#'   Default is \code{NA}, which falls back to \code{stats::approx}.
#'   See Details for more information.
#' @param xout_parameter_name character. The name of the parameter used
#'   by \code{interpolation_func} to receive the prediction vector.
#'   Default is \code{"xout"}. See Details for more information.
#' @param ... Additional parameters passed forward to
#'   \code{\link{ATL03_ATL08_photons_seg_dt_fitground}}.
#'
#' @return Returns the input \code{atl03_atl08_seg_dt} object with the
#'   \code{ph_h} column updated in place to contain height above ground
#'   level (AGL) in meters, computed as \code{h_ph - ground_elevation}.
#'   Values are \code{NA} for photons outside the interpolated ground range.
#'
#' @details
#' This function is a wrapper around
#' \code{\link{ATL03_ATL08_photons_seg_dt_fitground}}. It first estimates
#' the ground elevation at every photon location using the smoothing and
#' interpolation approach described in that function, then subtracts the
#' estimated ground elevation from the raw photon height (\code{h_ph})
#' to produce height above ground level.
#'
#' For full details on the smoothing window, smoothing function, and
#' interpolation function behaviour, see
#' \code{\link{ATL03_ATL08_photons_seg_dt_fitground}}.
#'
#' @examples
#' # Specifying the path to ATL03 and ATL08 files
#' atl03_path <- system.file("extdata", "atl03_clip.h5", package = "ICESat2VegR")
#' atl08_path <- system.file("extdata", "atl08_clip.h5", package = "ICESat2VegR")
#'
#' # Reading ATL03 and ATL08 data (h5 files)
#' atl03_h5 <- ATL03_read(atl03_path = atl03_path)
#' atl08_h5 <- ATL08_read(atl08_path = atl08_path)
#'
#' # Extracting ATL03 and ATL08 photons and heights
#' atl03_atl08_dt <- ATL03_ATL08_photons_attributes_dt_join(atl03_h5, atl08_h5)
#'
#' # Converting to seg_dt class (required input)
#' atl03_atl08_seg_dt <- ATL03_ATL08_segment_create(
#'   atl03_atl08_dt,
#'   segment_length = 30
#' )
#'
#' # Normalize photon heights using default ATBD smoothing window
#' atl03_atl08_seg_dt_norm <- ATL03_ATL08_photons_seg_dt_height_normalize(
#'   atl03_atl08_seg_dt,
#'   interpolation_func = approx,
#'   xout_parameter_name = "xout"
#' )
#' head(atl03_atl08_seg_dt_norm)
#'
#' # Normalize using a fixed 5 m smoothing window and mean aggregation
#' atl03_atl08_seg_dt_norm2 <- ATL03_ATL08_photons_seg_dt_height_normalize(
#'   atl03_atl08_seg_dt,
#'   smoothing_window = 5,
#'   smoothing_func = mean,
#'   interpolation_func = approx,
#'   xout_parameter_name = "xout"
#' )
#' head(atl03_atl08_seg_dt_norm2)
#'
#' close(atl03_h5)
#' close(atl08_h5)
#'
#' @import data.table
#' @export
ATL03_ATL08_photons_seg_dt_height_normalize <- function(
    atl03_atl08_seg_dt,
    smoothing_window = NA,
    smoothing_func = median,
    interpolation_func = NA,
    xout_parameter_name = "xout",
    ...) {
  ph_h <- dist_ph_along <- classed_pc_flag <- h_ph <- NA

  stopifnot(
    "atl03_atl08_seg_dt seems to be invalid, use the package function" =
      inherits(atl03_atl08_seg_dt, "icesat2.atl03_atl08_seg_dt")
  )

  if (!inherits(interpolation_func, "function")) {
    interpolation_func <- stats::approx
  }

  params <- list(
    ...
  )
  params[[xout_parameter_name]] <- list(atl03_atl08_seg_dt[, list(dist_ph_along)])


  elevation <- ATL03_ATL08_photons_seg_dt_fitground(
    atl03_atl08_seg_dt,
    smoothing_window,
    smoothing_func,
    interpolation_func,
    xout_parameter_name,
    ...
  )


  atl03_atl08_seg_dt[
    ,
    ph_h := h_ph - elevation
  ]

  range(atl03_atl08_seg_dt[classed_pc_flag >= 1, list(ph_h)], na.rm = TRUE)

  atl03_atl08_seg_dt
}
