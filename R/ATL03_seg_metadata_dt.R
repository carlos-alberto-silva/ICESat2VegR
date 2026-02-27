# Map of ATL03 geolocation (segment-level) variables
ATL03.seg.map <- list()
ATL03.seg.map[["altitude_sc"]]            <- "geolocation/altitude_sc"
ATL03.seg.map[["beta_angle"]]             <- "geolocation/beta_angle"
ATL03.seg.map[["bounce_time_offset"]]     <- "geolocation/bounce_time_offset"
ATL03.seg.map[["delta_time"]]             <- "geolocation/delta_time"
ATL03.seg.map[["full_sat_fract"]]         <- "geolocation/full_sat_fract"
ATL03.seg.map[["near_sat_fract"]]         <- "geolocation/near_sat_fract"
ATL03.seg.map[["neutat_delay_derivative"]]<- "geolocation/neutat_delay_derivative"
ATL03.seg.map[["neutat_delay_total"]]     <- "geolocation/neutat_delay_total"
ATL03.seg.map[["neutat_ht"]]              <- "geolocation/neutat_ht"
ATL03.seg.map[["ph_index_beg"]]           <- "geolocation/ph_index_beg"
ATL03.seg.map[["pitch"]]                  <- "geolocation/pitch"
ATL03.seg.map[["podppd_flag"]]            <- "geolocation/podppd_flag"
ATL03.seg.map[["range_bias_corr"]]        <- "geolocation/range_bias_corr"
ATL03.seg.map[["ref_azimuth"]]            <- "geolocation/ref_azimuth"
ATL03.seg.map[["ref_elev"]]               <- "geolocation/ref_elev"
ATL03.seg.map[["reference_photon_index"]] <- "geolocation/reference_photon_index"
ATL03.seg.map[["reference_photon_lat"]]   <- "geolocation/reference_photon_lat"
ATL03.seg.map[["reference_photon_lon"]]   <- "geolocation/reference_photon_lon"
ATL03.seg.map[["roll"]]                   <- "geolocation/roll"
ATL03.seg.map[["segment_dist_x"]]         <- "geolocation/segment_dist_x"
ATL03.seg.map[["segment_id"]]             <- "geolocation/segment_id"
ATL03.seg.map[["segment_length"]]         <- "geolocation/segment_length"
ATL03.seg.map[["segment_ph_cnt"]]         <- "geolocation/segment_ph_cnt"
ATL03.seg.map[["sigma_across"]]           <- "geolocation/sigma_across"
ATL03.seg.map[["sigma_along"]]            <- "geolocation/sigma_along"
ATL03.seg.map[["sigma_h"]]                <- "geolocation/sigma_h"
ATL03.seg.map[["sigma_lat"]]              <- "geolocation/sigma_lat"
ATL03.seg.map[["sigma_lon"]]              <- "geolocation/sigma_lon"
ATL03.seg.map[["solar_azimuth"]]          <- "geolocation/solar_azimuth"
ATL03.seg.map[["solar_elevation"]]        <- "geolocation/solar_elevation"
ATL03.seg.map[["surf_type"]]              <- "geolocation/surf_type"
ATL03.seg.map[["tx_pulse_energy"]]        <- "geolocation/tx_pulse_energy"
ATL03.seg.map[["tx_pulse_skew_est"]]      <- "geolocation/tx_pulse_skew_est"
ATL03.seg.map[["tx_pulse_width_lower"]]   <- "geolocation/tx_pulse_width_lower"
ATL03.seg.map[["tx_pulse_width_upper"]]   <- "geolocation/tx_pulse_width_upper"
ATL03.seg.map[["velocity_sc"]]            <- "geolocation/velocity_sc"
ATL03.seg.map[["yaw"]]                    <- "geolocation/yaw"

#' ATL03 geolocation segment metadata
#'
#' @description
#' Extract geolocation segment-level metadata from ICESat-2 ATL03 data.
#' Each row in the output corresponds to a single **20m geolocation segment**
#' along track, not to individual photons.
#'
#' In addition to variables from the ATL03 `geolocation` group, this
#' function derives a segment-level reference photon height (`h_ph`) using
#' the reference photon index and the photon-height array in the
#' `heights` group.
#'
#' @param atl03_h5
#'   An ICESat-2 ATL03 object (output of [ATL03_read()]),
#'   i.e. an S4 object of class
#'   [`ICESat2VegR::icesat2.atl03_h5`].
#' @param beam
#'   Character vector indicating beams to process
#'   (e.g. `gt1l`, `gt1r`, `gt2l`, `gt2r`, `gt3l`, `gt3r`).
#'   Only beams present in `atl03_h5$beams` are used.
#' @param attributes
#'   Character vector naming the segment-level variables to extract.
#'   By default, a broad set of geolocation and quality variables is used,
#'   including a derived reference-photon height (`h_ph`).
#'
#' @return
#' A [`data.table::data.table`] with one row per ATL03 geolocation segment,
#' with class `icesat2.atl03_seg_dt` prepended. Columns include:
#'
#' \itemize{
#'   \item `beam` - beam ID (e.g., `gt2l`).
#'   \item `strong_beam` - logical flag indicating whether the beam is
#'         classified as strong for the orbit.
#'   \item selected geolocation fields (see below).
#'   \item a derived `h_ph` column: height of the reference photon
#'         for each segment (one value per segment).
#' }
#'
#' @details
#' The following variables may be requested via `attributes``:
#'
#' - `h_ph`: Height of the **reference photon** above the WGS84 ellipsoid
#'   for each geolocation segment. This is derived from
#'   `geolocation/ph_index_beg`, `geolocation/reference_photon_index`,
#'   and `heights/h_ph`.
#' - `altitude_sc`: Height of the spacecraft above the WGS84 ellipsoid.
#' - `beta_angle`: Acute angle between Sun vector and orbit
#'    plane
#' - `bounce_time_offset`: Difference between the transmit time and the
#'   ground-bounce time of the reference photon.
#' - `delta_time`: Transmit time of the reference photon, measured in
#'   seconds from `atlas_sdp_gps_epoch``.
#' - `full_sat_fract`: Fraction of pulses within the segment that are
#'   fully saturated.
#' - `near_sat_fract`: Fraction of pulses within the segment that are
#'   nearly saturated.
#' - `neutat_delay_derivative`: Change in neutral atmospheric delay per
#'   unit height change.
#' - `neutat_delay_total`: Total neutral atmosphere delay correction
#'   (wet + dry).
#' - `neutat_ht`: Reference height of the neutral atmosphere range
#'   correction.
#' - `ph_index_beg`: 1-based index of the first photon in this segment
#'   within the photon-rate data.
#' - `pitch`: Spacecraft pitch (degrees), 3-2-1 Euler sequence.
#' - `podppd_flag`: Composite flag describing the quality of input
#'   geolocation products for the segment.
#' - `range_bias_corr`: Estimated range bias from geolocation analysis.
#' - `ref_azimuth`: Azimuth (radians) of the unit pointing vector for the
#'   reference photon in the local ENU frame.
#' - `ref_elev`: Elevation (radians) of the unit pointing vector for the
#'   reference photon in the local ENU frame.
#' - `reference_photon_index`: Index of the reference photon within the
#'   photon set for a segment.
#' - `reference_photon_lat`: Latitude of the reference photon.
#' - `reference_photon_lon`: Longitude of the reference photon.
#' - `roll`: Spacecraft roll (degrees), 3-2-1 Euler sequence.
#' - `segment_dist_x`: Along-track distance from the equator crossing to
#'   the start of the 20 m geolocation segment.
#' - `segment_id`: 7-digit along-track geolocation segment identifier.
#' - `segment_length`: Along-track length of the geolocation segment
#'   (typically 20 m).
#' - `segment_ph_cnt`: Number of photons in the segment.
#' - `sigma_across`: Estimated Cartesian across-track uncertainty
#'   (1-sigma) for the reference photon.
#' - `sigma_along`: Estimated Cartesian along-track uncertainty (1-sigma)
#'   for the reference photon.
#' - `sigma_h`: Estimated height uncertainty (1-sigma) for the reference
#'   photon bounce point.
#' - `sigma_lat`: Estimated geodetic latitude uncertainty (1-sigma) for
#'   the reference photon.
#' - `sigma_lon`: Estimated geodetic longitude uncertainty (1-sigma) for
#'   the reference photon.
#' - `solar_azimuth`: Azimuth (degrees east) of the sun position vector
#'   from the reference photon bounce point in the local ENU frame.
#' - `solar_elevation`: Elevation (degrees) of the sun position vector
#'   from the reference photon bounce point in the local ENU frame.
#' - `surf_type`: Flags describing which surface types the segment is
#'   associated with (land, ocean, sea ice, land ice, inland water).
#' - `tx_pulse_energy`: Average transmit pulse energy per beam.
#' - `tx_pulse_skew_est`: Difference between the averages of the lower and
#'   upper threshold crossing times, estimating transmit pulse skew.
#' - `tx_pulse_width_lower`: Average distance between lower threshold
#'   crossing times measured by the Start Pulse Detector.
#' - `tx_pulse_width_upper`: Average distance between upper threshold
#'   crossing times measured by the Start Pulse Detector.
#' - `velocity_sc`: Spacecraft velocity components (east
#'   component, north component, up component)
#'   an observer on the ground would measure.
#'   While values are common to all beams, this
#'   parameter is naturally produced as part of
#'   geolocation.
#' - `yaw`: Spacecraft yaw (degrees), 3-2-1 Euler sequence.
#'
#' @seealso
#' \url{https://icesat-2.gsfc.nasa.gov/sites/default/files/page_files/ICESat2_ATL03_ATBD_r006.pdf}
#'
#' @examples
#' \dontrun{
#' atl03_path <- system.file(
#'   "extdata", "atl03_clip.h5",
#'   package = "ICESat2VegR"
#' )
#'
#' atl03_h5 <- ATL03_read(atl03_path)
#'
#' # Extract ATL03 geolocation segment metadata
#' atl03_segment_dt <- ATL03_seg_metadata_dt(atl03_h5)
#'
#' head(atl03_segment_dt)
#' close(atl03_h5)
#' }
#'
#' @export
ATL03_seg_metadata_dt <- function(
  atl03_h5,
  beam = c("gt1l", "gt1r", "gt2l", "gt2r", "gt3l", "gt3r"),
  attributes = c(
    "h_ph",
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
  )
) {
  # Input check
  if (!inherits(atl03_h5, "icesat2.atl03_h5")) {
    stop("atl03_h5 must be an object of class 'icesat2.atl03_h5' (output of ATL03_read()).")
  }

  `:=` <- data.table::`:=`
  strong_beam <- h_ph <- NA  # for R CMD check

  # Restrict beams to those present
  beam <- intersect(atl03_h5$beams, beam)

  seg.dt <- list()
  pb <- utils::txtProgressBar(min = 0, max = length(beam), style = 3, file = stderr())
  i_s <- 0

  # We always need reference_photon_lon / lat for surf_type masking logic
  mask_surf_type <- attributes == "surf_type"
  attributes <- c(attributes[!mask_surf_type], c("reference_photon_lon", "reference_photon_lat"))

  # Loop over beams
  for (ii in beam) {
    i_s <- i_s + 1
    beam_group <- atl03_h5[[ii]]
    dt <- data.table::data.table()

    # Loop over requested attributes
    for (attr in attributes) {
      dataset_name <- ATL03.seg.map[[attr]]

      if (!is.null(dataset_name) && beam_group$exists(dataset_name)) {
        # 1D datasets
        if (length(beam_group[[dataset_name]]$dims) == 1) {
          dt[, (attr) := beam_group[[dataset_name]][]]
        # 2D datasets: expand across first dimension
        } else if (length(beam_group[[dataset_name]]$dims) == 2) {
          for (jj in seq_len(beam_group[[dataset_name]]$dims[1])) {
            attr_num <- gettextf("%s_%s", attr, jj)
            dt[, (attr_num) := beam_group[[dataset_name]][jj, ]]
          }
        }
      } else {
        # fall back: direct name under beam group
        if (beam_group$exists(attr)) {
          dt[, (attr) := beam_group[[attr]][]]
        }
      }

      # Derived reference-photon height at segment scale
      if (attr == "h_ph") {
        ref_idx <- beam_group[["geolocation/reference_photon_index"]][]
        ref_idx_mask <- ref_idx > 0
        idx_mask <- seq_along(ref_idx_mask)[!ref_idx_mask]

        ph_index_beg <- beam_group[["geolocation/ph_index_beg"]][]

        if (beam_group[["heights/h_ph"]]$dims > 0) {
          reference_idx <- ph_index_beg + ref_idx - 1
          reference_idx[idx_mask] <- 1

          dt[, h_ph := beam_group[["heights/h_ph"]][reference_idx]]
          dt[idx_mask, h_ph := NA_real_]
        }
      }
    }

    dt[, beam := ii]

    if (ncol(dt) > 1) {
      seg.dt[[i_s]] <- dt
    }

    utils::setTxtProgressBar(pb, i_s)
  }

  seg.dt <- data.table::rbindlist(seg.dt, fill = TRUE)
  seg.dt[["strong_beam"]] <- seg.dt$beam %in% atl03_h5$strong_beams
  seg.dt <- stats::na.omit(seg.dt)

  prepend_class(seg.dt, "icesat2.atl03_seg_dt")

  close(pb)
  seg.dt
}
