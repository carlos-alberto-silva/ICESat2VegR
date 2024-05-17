#' Fit and estimate ground elevation for photons or arbitrary distances
#' from the track beginning.
#'
#' @description
#' Function to estimate ground elevation using smoothing and
#' interpolation functions
#'
#' @param atl03_atl08_seg_dt An S4 object of class [`ICESat2VegR::icesat2.atl03_atl08_seg_dt-class`] containing ATL03 and ATL08 data
#' (output of [ATL03_ATL08_photons_attributes_dt_join()] function).
#' @param smoothing_window numeric. The smoothing window size in meters for smoothing the photon cloud.
#' Default is NA, see details for more information.
#' @param smoothing_func function. The smoothing function to be applied on the smoothing window.
#' @param interpolation_func function. The interpolation function to estimate the ground elevation.
#' Default [`stats::approx()`].
#' @param xout_parameter_name character. The parameter name used for the `interpolation_func` to
#' which will be used to predict, default "xout".
#' @param ... parameters to be passed forward to [`ATL03_ATL08_photons_seg_dt_fitground()`].
#'
#' @details
#' The function for calculating the ground will first pass a smoothing
#' window with `smoothing_window` size,
#' applying the `smoothing_func` to aggregate the ground photons.
#'
#' Then it will use an interpolation function between those
#' aggregated photons to calculate a smooth surface.
#'
#' The `smoothing_func` signature will depend on the function used.
#' It is assumed that the first two arguments are vectors of `x` (independent
#' variable) and `y` (the prediction to be interpolated). The remaining
#' arguments are passed through `...`.
#'
#' The interpolation functions need a third parameter which is the
#' `x` vector to be interpolated. Functions from `stats` base package
#' `stats::approx()` and `stats::spline()` name this argument as `xout`,
#' so you can use:
#'
#' ```
#' ATL03_ATL08_photons_fitground_seg_dt(
#'   dt,
#'   interpolation_func = approx,
#'   xout = 1:30
#' )
#' ```
#'
#' For example, to interpolate the values for the 1:30 vector. But other functions
#' may name the parameter differently, such as [`signal::pchip()`], which name
#' the parameter as `xi` instead of `xout`. [`signal::pchip()`] is the
#' algorithm used by ATL08 ATBD.
#'
#' The `smoothing_window` can be left NA, which will use the ATBD algoritm
#' for calculating the window size:
#'
#' \eqn{Sspan = ceil[5 + 46 * (1 - e^{-a * length})]}, where *length*
#' is the number of photons within segment.
#'
#' \deqn{a \approx 21x10^{-6}}
#'
#' \deqn{window_size = \frac{2}{3} Sspan}
#'
#' This is not the same algorithm as used in ATL08
#' but is an adapted version that uses the ATL08
#' pre-classification
#' @export
ATL03_ATL08_photons_dt_height_normalize <- function(
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

  range(atl03_atl08_seg_dt[classed_pc_flag >= 1, list(ph_h)], na.rm = T)

  atl03_atl08_seg_dt
}
