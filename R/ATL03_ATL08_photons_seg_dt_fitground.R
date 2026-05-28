#' Fit and estimate ground elevation for photons or arbitrary distances
#' from the track beginning
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
#' @param xout_parameter_name character. Optional, can be used to inform the parameter name that the
#' interpolation_func uses for passing the prediction vector and already use the photons for prediction.
#' Default NA will use the ...
#' @param ... Optional parameters to pass to the interpolation_func, see details for more information.
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
#' For example, to interpolate the values for the 1:30 vector. However, other
#' functions may name the parameter differently, such as `signal::pchip()`,
#' which calls the parameter `xi` instead of `xout`.
#' The `pchip` algorithm (as implemented in the \strong{signal} package)
#' is the one used by the ATL08 ATBD.
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
ATL03_ATL08_photons_seg_dt_fitground <- function(
    atl03_atl08_seg_dt,
    smoothing_window = NA,
    smoothing_func = median,
    interpolation_func = NA,
    xout_parameter_name = "xout",
    ...) {
  classed_pc_flag <-
    dist_ph_along <-
    segment_id <-
    h_ph <- NA

  stopifnot(
    "atl03_atl08_seg_dt seems to be invalid, use the package function" =
      inherits(atl03_atl08_seg_dt, "icesat2.atl03_atl08_seg_dt")
  )

  if (!inherits(interpolation_func, "function")) {
    interpolation_func <- approxfun
  }

  windowSize <- NA
  if (!is.numeric(smoothing_window)) {
    windowSize <- atl03_atl08_seg_dt[
      classed_pc_flag > 0,
      list(smoothing_window = 2 / 3 * ceiling(5 + 46 * (1 - exp(-21e-6 * .N)))),
      by = segment_id
    ]
  }

  ground_photons <- atl03_atl08_seg_dt[
    classed_pc_flag == 1 # Use only ground photons
  ]


  if (!is.numeric(smoothing_window)) {
    ground_photons <- ground_photons[windowSize, , on = "segment_id"]
  }

  smoothed <- ground_photons[
    ,
    list(h_ph = smoothing_func(h_ph)),
    by = list(
      dist_ph_along = floor(
        dist_ph_along / smoothing_window
      ) * smoothing_window + (smoothing_window / 2)
    ),
  ][order(dist_ph_along)][
    !is.na(dist_ph_along),
    list(
      dist_ph_along,
      h_ph
    )
  ]

  if (is.na(xout_parameter_name)) {
    return(interpolation_func(
      smoothed$dist_ph_along,
      smoothed$h_ph,
      ...
    ))
  } else {
    params <- list(
      smoothed$dist_ph_along,
      smoothed$h_ph,
      ...
    )
    params[[xout_parameter_name]] <-
      atl03_atl08_seg_dt$dist_ph_along

    return(do.call(
      interpolation_func,
      params
    ))
  }
}
