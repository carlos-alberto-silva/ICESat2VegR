# ATL08.var.map
ATL03.seg.map <- list()
ATL03.seg.map[["altitude_sc"]] <- "geolocation/altitude_sc"
ATL03.seg.map[["bounce_time_offset"]] <- "geolocation/bounce_time_offset"
ATL03.seg.map[["delta_time"]] <- "geolocation/delta_time"
ATL03.seg.map[["full_sat_fract"]] <- "geolocation/full_sat_fract"
ATL03.seg.map[["near_sat_fract"]] <- "geolocation/near_sat_fract"
ATL03.seg.map[["neutat_delay_derivative"]] <- "geolocation/neutat_delay_derivative"
ATL03.seg.map[["neutat_delay_total"]] <- "geolocation/neutat_delay_total"
ATL03.seg.map[["neutat_ht"]] <- "geolocation/neutat_ht"
ATL03.seg.map[["ph_index_beg"]] <- "geolocation/ph_index_beg"
ATL03.seg.map[["pitch"]] <- "geolocation/pitch"
ATL03.seg.map[["podppd_flag"]] <- "geolocation/podppd_flag"
ATL03.seg.map[["range_bias_corr"]] <- "geolocation/range_bias_corr"
ATL03.seg.map[["ref_azimuth"]] <- "geolocation/ref_azimuth"
ATL03.seg.map[["ref_elev"]] <- "geolocation/ref_elev"
ATL03.seg.map[["reference_photon_index"]] <- "geolocation/reference_photon_index"
ATL03.seg.map[["reference_photon_lat"]] <- "geolocation/reference_photon_lat"
ATL03.seg.map[["reference_photon_lon"]] <- "geolocation/reference_photon_lon"
ATL03.seg.map[["roll"]] <- "geolocation/roll"
ATL03.seg.map[["segment_dist_x"]] <- "geolocation/segment_dist_x"
ATL03.seg.map[["segment_id"]] <- "geolocation/segment_id"
ATL03.seg.map[["segment_length"]] <- "geolocation/segment_length"
ATL03.seg.map[["segment_ph_cnt"]] <- "geolocation/segment_ph_cnt"
ATL03.seg.map[["sigma_across"]] <- "geolocation/sigma_across"
ATL03.seg.map[["sigma_along"]] <- "geolocation/sigma_along"
ATL03.seg.map[["sigma_h"]] <- "geolocation/sigma_h"
ATL03.seg.map[["sigma_lat"]] <- "geolocation/sigma_lat"
ATL03.seg.map[["sigma_lon"]] <- "geolocation/sigma_lon"
ATL03.seg.map[["solar_azimuth"]] <- "geolocation/solar_azimuth"
ATL03.seg.map[["solar_elevation"]] <- "geolocation/solar_elevation"
ATL03.seg.map[["surf_type"]] <- "geolocation/surf_type"
ATL03.seg.map[["tx_pulse_energy"]] <- "geolocation/tx_pulse_energy"
ATL03.seg.map[["tx_pulse_skew_est"]] <- "geolocation/tx_pulse_skew_est"
ATL03.seg.map[["tx_pulse_width_lower"]] <- "geolocation/tx_pulse_width_lower"
ATL03.seg.map[["tx_pulse_width_upper"]] <- "geolocation/tx_pulse_width_upper"
ATL03.seg.map[["yaw"]] <- "geolocation/yaw"

#' ATL03 segments attributes
#'
#' @description This function extracts segment attributes from ICESat-2 ATL03 data
#'
#' @param atl03_h5 A ICESat-2 ATL03 object (output of [ATL03_read()] function).
#' An S4 object of class [`ICESat2VegR::icesat2.atl03_dt-class`].
#' @param beam Character vector indicating beams to process (e.g. "gt1l", "gt1r", "gt2l", "gt2r", "gt3l", "gt3r")
#' @param attributes Character vector indicating the attrivutes
#'
#' @return Returns an S4 object of class [data.table::data.table]
#' containing the ATL03 segment attributes.
#'
#' @details These are the available variables for extraction:
#'
#' - `altitude_sc`: Height of the spacecraft above the WGS84 ellipsoid.
#' - `bounce_time_offset`: The difference between the transmit time and the ground bounce time of the reference photons.
#' - `delta_time`: Transmit time of the reference photon, measured in seconds from the atlas_sdp_gps_epoch.
#' - `full_sat_fract`: The fraction of pulses within the segment determined to be fully saturated.
#' - `near_sat_fract`: The fraction of pulses within the segment determined to be nearly saturated.
#' - `neutat_delay_derivative`: Change in neutral atmospheric delay per height change.
#' - `neutat_delay_total`: Total neutral atmosphere delay correction (wet+dry).
#' - `neutat_ht`: Reference height of the neutral atmosphere range correction.
#' - `ph_index_beg`: Index (1-based) within the photon-rate data of the first photon within this segment.
#' - `pitch`: Spacecraft pitch, computed using a 3, 2, 1 Euler angle sequence, with units in degrees.
#' - `podppd_flag`: A composite flag indicating the quality of input geolocation products for the specific ATL03 segment.
#' - `range_bias_corr`: The estimated range bias from geolocation analysis.
#' - `ref_azimuth`: Azimuth of the unit pointing vector for the reference photon in the local ENU frame, in radians.
#' - `ref_elev`: Elevation of the unit pointing vector for the reference photon in the local ENU frame, in radians.
#' - `reference_photon_index`: Index of the reference photon within the set of photons grouped within a segment.
#' - `reference_photon_lat`: Latitude of each reference photon.
#' - `reference_photon_lon`: Longitude of each reference photon.
#' - `roll`: Spacecraft roll, computed using a 3, 2, 1 Euler angle sequence, with units in degrees.
#' - `segment_dist_x`: Along-track distance from the equator crossing to the start of the 20 meter geolocation segment.
#' - `segment_id`: A 7-digit number identifying the along-track geolocation segment number.
#' - `segment_length`: The along-track length of the along-track segment, typically 20 meters.
#' - `segment_ph_cnt`: Number of photons in a given along-track segment.
#' - `sigma_across`: Estimated Cartesian across-track uncertainty (1-sigma) for the reference photon.
#' - `sigma_along`: Estimated cartesian along-track uncertainty (1-sigma) for the reference photon.
#' - `sigma_h`: Estimated height uncertainty (1-sigma) for the reference photon bounce point.
#' - `sigma_lat`: Estimated geodetic Latitude uncertainty (1-sigma) for the reference photon bounce point.
#' - `sigma_lon`: Estimated geodetic Longitude uncertainty (1-sigma) for the reference photon bounce point.
#' - `solar_azimuth`: Azimuth of the sun position vector from the reference photon bounce point position in
#' the local ENU frame, in degrees east.
#' - `solar_elevation`: Elevation of the sun position vector from the reference photon bounce point
#' position in the local ENU frame, in degrees.
#' - `surf_type`: Flags describing which surface types an interval is associated with (e.g., land,
#' ocean, sea ice, land ice, inland water).
#' - `tx_pulse_energy`: Average transmit pulse energy, measured by the internal laser energy monitor,
#' split into per-beam measurements.
#' - `tx_pulse_skew_est`: Difference between the averages of the lower and upper threshold crossing times,
#' estimating the transmit pulse skew.
#' - `tx_pulse_width_lower`: Average distance between the lower threshold crossing times measured by the
#' Start Pulse Detector.
#' - `tx_pulse_width_upper`: Average distance between the upper threshold crossing times measured by the
#' Start Pulse Detector.
#' - `velocity_sc`: Spacecraft velocity components (east component, north component, up component)
#' an observer on the ground would measure.
#' - `yaw`: Spacecraft yaw, computed using a 3, 2, 1 Euler angle sequence, with units in degrees.
#'
#' @seealso \url{https://icesat-2.gsfc.nasa.gov/sites/default/files/page_files/ICESat2_ATL03_ATBD_r006.pdf}
#'
#' @examples
#' # Specifying the path to ATL03 file
#' atl03_path <- system.file("extdata",
#'   "atl03_clip.h5",
#'   package = "ICESat2VegR"
#' )
#'
#' # Reading ATL03 data (h5 file)
#' atl03_h5 <- ATL03_read(atl03_path = atl03_path)
#'
#' # Extracting ATL03 segment attributes
#' atl03_segment_dt <- ATL03_seg_attributes_dt(atl03_h5 = atl03_h5)
#'
#' head(atl03_segment_dt)
#' close(atl03_h5)
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
  if (!inherits(atl03_h5, "icesat2.atl03_h5")) {
    stop("atl03_h5 must be an object of class 'icesat2.atl03_h5' - output of [ATL03_read()] function ")
  }

  `:=` <- data.table::`:=`

  # Initialize data.table inside variable to avoid R CMD check warnings
  strong_beam <- h_ph <- NA

  # Check beams to select
  beam <- intersect(atl03_h5$beams, beam)

  seg.dt <- list()

  pb <- utils::txtProgressBar(min = 0, max = length(beam), style = 3)

  i_s <- 0

  mask_surf_type <- attributes == "surf_type"
  attributes <- c(attributes[!mask_surf_type], c("reference_photon_lon", "reference_photon_lat"))

  # Loop over each beam
  for (ii in beam) {
    i_s <- i_s + 1 # Increment counter
    beam_group <- atl03_h5[[ii]] # Extract ith beam group from ATL03 data
    dt <- data.table::data.table() # Initialize data table

    # Loop within each beam to process attributes
    for (attr in attributes) {
      dataset_name <- ATL03.seg.map[[attr]] # Get dataset name for attribute

      # If dataset exists, handle different dimensions
      if (!is.null(dataset_name) && beam_group$exists(dataset_name)) {
        # For dimension=1, assign directly
        if (length(beam_group[[dataset_name]]$dims) == 1) {
          dt[, eval(attr) := beam_group[[dataset_name]][]]
        }
        # For dimension=2, create extra columns per dimension
        else if (length(beam_group[[dataset_name]]$dims) == 2) {
          for (jj in seq_len(beam_group[[dataset_name]]$dims[1])) {
            attr_num <- gettextf("%s_%s", attr, jj)
            dt[, eval(attr_num) := beam_group[[dataset_name]][jj, ]]
          }
        }
      } else {
        if (beam_group$exists(attr)) {
          dt[, eval(attr) := beam_group[[attr]][]]
        }
      }

      # Special case when attribute is "h_ph"
      if (attr == "h_ph") {
        ref_idx <- beam_group[["geolocation/reference_photon_index"]][] # Get reference photon index
        idx_mask <- seq_along(ref_idx)[ref_idx > 0] # Get valid indices (greater than 0)
        dt[idx_mask, h_ph := beam_group[["heights/h_ph"]][ref_idx[idx_mask]]] # Update h_ph for those indices
      }
    }

    dt[, beam := ii] # Add beam number as a new column to the data table

    # If there are more than one columns in data table, store it into list 'seg.dt'
    if (ncol(dt) > 1) {
      seg.dt[[i_s]] <- dt
    }

    utils::setTxtProgressBar(pb, i_s) # Update progress bar
  }

  seg.dt <- data.table::rbindlist(seg.dt, fill = TRUE)
  seg.dt[, strong_beam := ifelse(beam %in% atl03_h5$strong_beams, TRUE, FALSE)]
  seg.dt <- na.omit(seg.dt)
  prepend_class(seg.dt, "icesat2.atl03_seg_dt")


  close(pb)

  return(seg.dt)
}
