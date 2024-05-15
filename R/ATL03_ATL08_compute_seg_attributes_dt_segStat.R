#' Statistics of ATL03 and ATL08 labeled photons at the segment level
#'
#' @description Computes a series of statistics from ATL03 and ATL08 labeled photons
#' within a given segment length.
#'
#' @param atl03_atl08_seg_dt  An S4 object of class [`ICESat2VegR::icesat2.atl03_atl08_seg_dt-class`] containing ATL03 and ATL08 data
#' (output of [ATL03_ATL08_photons_attributes_dt_join()] function).
#' @param list_expr The function to be applied for computing the defined statistics
#' @param seg_length Segment length. Default is 30 m
#' @param ph_class Character vector indicating photons to process based
#' on the classification (1=ground, 2=canopy, 3=top canopy),
#' Default is c(2,3)
#' @param beam Character vector indicating beams to process. Default is
#' c("gt1l", "gt1r", "gt2l", "gt2r", "gt3l", "gt3r")
#' @param quality_ph Indicates the quality of the associated photon.
#' 0 = nominal, 1 = possible_afterpulse, 2 = possible_impulse_response_
#' effect, 3=possible_tep. Default is 0
#' @param night_flag Flag indicating the data were acquired in night conditions: 0=day, 1=night. Default is 1
#'
#' @return Returns an S4 object of class [`ICESat2VegR::icesat2.atl08_dt-class`]
#' Containing Statistics of ATL03 and ATL08 labeled photons
#'
#' @examples
#' # Specifying ATL03 and ATL08 file path
#' atl03_path <- system.file("extdata",
#'   "atl03_clip.h5",
#'   package = "ICESat2VegR"
#' )
#'
#' atl08_path <- system.file("extdata",
#'   "atl08_clip.h5",
#'   package = "ICESat2VegR"
#' )
#'
#' # Reading ATL03 data (h5 file)
#' atl03_h5 <- ATL03_read(atl03_path = atl03_path)
#'
#' # Reading ATL08 data (h5 file)
#' atl08_h5 <- ATL08_read(atl08_path = atl08_path)
#'
#' # Extracting ATL03 and ATL08 labeled photons
#' atl03_atl08_dt <- ATL03_ATL08_photons_attributes_dt_join(atl03_h5, atl08_h5)
#'
#' # Computing the max canopy height at 30 m segments
#' atl03_atl08_dt_seg <- ATL03_ATL08_segment_create(atl03_atl08_dt, segment_length = 30)
#'
#' max_canopy <- ATL03_ATL08_compute_seg_attributes_dt_segStat(atl03_atl08_dt_seg,
#'   list_expr = max(ph_h),
#'   ph_class = c(2, 3),
#'   beam = c("gt1l", "gt1r", "gt2l", "gt2r", "gt3l", "gt3r"),
#'   quality_ph = 0,
#'   night_flag = 0
#' )
#'
#' head(max_canopy)
#'
#' # Computing a series of canopy height statistics from customized list expressions
#' canopy_metrics <- ATL03_ATL08_compute_seg_attributes_dt_segStat(atl03_atl08_dt_seg,
#'   list_expr = list(
#'     max_ph_elevation = max(h_ph),
#'     h_canopy = quantile(ph_h, 0.98),
#'     n_canopy = sum(classed_pc_flag == 2),
#'     n_top_canopy = sum(classed_pc_flag == 3)
#'   ),
#'   ph_class = c(2, 3),
#'   beam = c("gt1l", "gt1r", "gt2l", "gt2r", "gt3l", "gt3r"),
#'   quality_ph = 0,
#'   night_flag = 0 # there are no night photons in this dataset
#' )
#'
#' head(canopy_metrics)
#'
#' close(atl03_h5)
#' close(atl08_h5)
#' @include lazy_applier.R
#' @import data.table lazyeval
#' @export
ATL03_ATL08_compute_seg_attributes_dt_segStat <- function(
    atl03_atl08_seg_dt,
    list_expr,
    ph_class = c(0, 1, 2, 3),
    beam = c("gt1l", "gt1r", "gt2l", "gt2r", "gt3l", "gt3r"),
    quality_ph = NULL,
    night_flag = NULL) {
  if (!inherits(atl03_atl08_seg_dt, "icesat2.atl03_atl08_seg_dt")) {
    stop("atl03_atl08_dt needs to be an object of class 'icesat2.atl03_atl08_seg_dt' ")
  }

  selected_quality_ph <- quality_ph
  selected_night_flag <- night_flag
  selected_beams <- beam
  classed_pc_flag <- NA

  atl03_atl08_seg_dt2 <- atl03_atl08_seg_dt[
    classed_pc_flag %in% ph_class &
      ifelse(is.null(selected_quality_ph), TRUE, quality_ph == selected_quality_ph) &
      beam %in% selected_beams &
      ifelse(is.null(selected_night_flag), TRUE, night_flag == selected_night_flag)
  ]

  args <- substitute(list_expr)

  metrics <- lazy_apply_dt_call(
    atl03_atl08_seg_dt2,
    args,
    "by = .(segment_id, beam, longitude = centroid_x, latitude = centroid_y)"
  )

  prepend_class(metrics, "icesat2.atl03_atl08_seg_dt")

  return(metrics)
}
