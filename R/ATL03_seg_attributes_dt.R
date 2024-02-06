# ATL08.var.map
ATL03.seg.map <- list()
ATL03.seg.map[["altitude_sc"]] <- "/geolocation/altitude_sc"
ATL03.seg.map[["bounce_time_offset"]] <- "/geolocation/bounce_time_offset"
ATL03.seg.map[["delta_time"]] <- "/geolocation/delta_time"
ATL03.seg.map[["full_sat_fract"]] <- "/geolocation/full_sat_fract"
ATL03.seg.map[["near_sat_fract"]] <- "/geolocation/near_sat_fract"
ATL03.seg.map[["neutat_delay_derivative"]] <- "/geolocation/neutat_delay_derivative"
ATL03.seg.map[["neutat_delay_total"]] <- "/geolocation/neutat_delay_total"
ATL03.seg.map[["neutat_ht"]] <- "/geolocation/neutat_ht"
ATL03.seg.map[["ph_index_beg"]] <- "/geolocation/ph_index_beg"
ATL03.seg.map[["pitch"]] <- "/geolocation/pitch"
ATL03.seg.map[["podppd_flag"]] <- "/geolocation/podppd_flag"
ATL03.seg.map[["range_bias_corr"]] <- "/geolocation/range_bias_corr"
ATL03.seg.map[["ref_azimuth"]] <- "/geolocation/ref_azimuth"
ATL03.seg.map[["ref_elev"]] <- "/geolocation/ref_elev"
ATL03.seg.map[["reference_photon_index"]] <- "/geolocation/reference_photon_index"
ATL03.seg.map[["reference_photon_lat"]] <- "/geolocation/reference_photon_lat"
ATL03.seg.map[["reference_photon_lon"]] <- "/geolocation/reference_photon_lon"
ATL03.seg.map[["roll"]] <- "/geolocation/roll"
ATL03.seg.map[["segment_dist_x"]] <- "/geolocation/segment_dist_x"
ATL03.seg.map[["segment_id"]] <- "/geolocation/segment_id"
ATL03.seg.map[["segment_length"]] <- "/geolocation/segment_length"
ATL03.seg.map[["segment_ph_cnt"]] <- "/geolocation/segment_ph_cnt"
ATL03.seg.map[["sigma_across"]] <- "/geolocation/sigma_across"
ATL03.seg.map[["sigma_along"]] <- "/geolocation/sigma_along"
ATL03.seg.map[["sigma_h"]] <- "/geolocation/sigma_h"
ATL03.seg.map[["sigma_lat"]] <- "/geolocation/sigma_lat"
ATL03.seg.map[["sigma_lon"]] <- "/geolocation/sigma_lon"
ATL03.seg.map[["solar_azimuth"]] <- "/geolocation/solar_azimuth"
ATL03.seg.map[["solar_elevation"]] <- "/geolocation/solar_elevation"
ATL03.seg.map[["surf_type"]] <- "/geolocation/surf_type"
ATL03.seg.map[["tx_pulse_energy"]] <- "/geolocation/tx_pulse_energy"
ATL03.seg.map[["tx_pulse_skew_est"]] <- "/geolocation/tx_pulse_skew_est"
ATL03.seg.map[["tx_pulse_width_lower"]] <- "/geolocation/tx_pulse_width_lower"
ATL03.seg.map[["tx_pulse_width_upper"]] <- "/geolocation/tx_pulse_width_upper"
ATL03.seg.map[["yaw"]] <- "/geolocation/yaw"

#' ATL03 segments attributes
#'
#' @description This function extracts segment attributes from ICESat-2 ATL03 data
#'
#' @param ATL03_h5 A ICESat-2 ATL03 object (output of [ATL03_read()] function).
#' An S4 object of class [ICESat2VegR::icesat2.atl03_dt].
#' @param beam Character vector indicating beams to process (e.g. "gt1l", "gt1r", "gt2l", "gt2r", "gt3l", "gt3r")
#' @param attributes Character vector indicating the attrivutes
#'
#' @return Returns an S4 object of class [data.table::data.table]
#' containing the ATL03 segment attributes.
#'
#' @details These are the available variables for extraction:
#' - `altitude_sc`
#' - `bounce_time_offset`
#' - `delta_time`
#' - `full_sat_fract`
#' - `near_sat_fract`
#' - `neutat_delay_derivative`
#' - `neutat_delay_total`
#' - `neutat_ht`
#' - `ph_index_beg`
#' - `pitch`
#' - `podppd_flag`
#' - `range_bias_corr`
#' - `ref_azimuth`
#' - `ref_elev`
#' - `reference_photon_index`
#' - `reference_photon_lat`
#' - `reference_photon_lon`
#' - `roll`
#' - `segment_dist_x`
#' - `segment_id`
#' - `segment_length`
#' - `segment_ph_cnt`
#' - `sigma_across`
#' - `sigma_along`
#' - `sigma_h`
#' - `sigma_lat`
#' - `sigma_lon`
#' - `solar_azimuth`
#' - `solar_elevation`
#' - `surf_type`
#' - `tx_pulse_energy`
#' - `tx_pulse_skew_est`
#' - `tx_pulse_width_lower`
#' - `tx_pulse_width_upper`
#' - `velocity_sc`
#' - `yaw`
#'
#' @seealso \url{https://icesat-2.gsfc.nasa.gov/sites/default/files/page_files/ICESat2_ATL03_ATBD_r006.pdf}
#'
#' @examples
#'
#' # Specifying the path to ATL03 file (zip file)
#' outdir <- tempdir()
#' atl03_zip <- system.file("extdata",
#'   "ATL03_20220401221822_01501506_005_01.zip",
#'   package = "ICESat2VegR"
#' )
#'
#' # Unzipping ATL03 file
#' atl03_path <- unzip(atl03_zip, exdir = outdir)
#'
#' # Reading ATL03 data (h5 file)
#' atl03_h5 <- ATL03_read(atl03_path = atl03_path)
#'
#' # Extracting ATL03 photons attributes
#' atl03_photons_dt <- ATL03_seg_attributes_dt(atl03_h5 = atl03_h5)
#' head(ATL03_photons_dt)
#'
#' close(ATL03_h5)
#' @export
ATL03_seg_attributes_dt <- function(atl03_h5,
                                        beam = c("gt1l", "gt1r", "gt2l", "gt2r", "gt3l", "gt3r"),
                                        attributes = c(
                                          "altitude_sc",
                                          "bounce_time_offset",
                                          "delta_time",
                                          "full_sat_fract",
                                          "near_sat_fract",
                                          "neutat_delay_derivative",
                                          "neutat_delay_total",
                                          "neutat_ht",
                                          "ph_index_beg",
                                          "pitch",
                                          "podppd_flag",
                                          "range_bias_corr",
                                          "ref_azimuth",
                                          "ref_elev",
                                          "reference_photon_index",
                                          "reference_photon_lat",
                                          "reference_photon_lon",
                                          "roll",
                                          "segment_dist_x",
                                          "segment_id",
                                          "segment_length",
                                          "segment_ph_cnt",
                                          "sigma_across",
                                          "sigma_along",
                                          "sigma_h",
                                          "sigma_lat",
                                          "sigma_lon",
                                          "solar_azimuth",
                                          "solar_elevation",
                                          "surf_type",
                                          "tx_pulse_energy",
                                          "tx_pulse_skew_est",
                                          "tx_pulse_width_lower",
                                          "tx_pulse_width_upper",
                                          "yaw"
                                        )) {
  # Check file input
  if (!class(atl03_h5) == "icesat2.atl03_h5") {
    stop("atl03_h5 must be an object of class 'icesat2.atl03_h5' - output of [ATL03_read()] function ")
  }

  # h5
  atl03_h5v2 <- atl03_h5@h5

  # Check beams to select
  groups_id <- getBeams(atl03_h5)

  check_beams <- groups_id %in% beam
  beam <- groups_id[check_beams]

  seg.dt <- list()

  pb <- utils::txtProgressBar(min = 0, max = length(beam), style = 3)

  i_s <- 0

  mask_surf_type <- attributes == "surf_type"
  has_surf_type = any(mask_surf_type)
  attributes <- attributes[! mask_surf_type]
  for (i in beam) {
    i_s <- i_s + 1

    
    dt <- data.table::data.table(sapply(attributes, function(x) {
      dataset_name <- sprintf('%s/%s', i, ATL03.seg.map[[x]])
      atl03_h5[[dataset_name]][]
    }))
    
    seg.dt[[i_s]] <- dt
    utils::setTxtProgressBar(pb, i_s)
  }

  seg.dt <- data.table::rbindlist(seg.dt)
  prepend_class(seg.dt, "icesat2.atl03_seg_dt")


  close(pb)

  return(seg.dt)
}
