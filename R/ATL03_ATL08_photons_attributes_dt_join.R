#' Join ATL03 and ATL08 photons attributes
#'
#' @description This function joins ATL03 and ATL08 computed photons attributes
#'
#' @param atl03_h5 A ICESat-2 ATL03 object (output of [ATL03_read()] function).
#' An S4 object of class [`ICESat2VegR::icesat2.atl03_dt-class`].
#' @param atl08_h5 A ICESat-2 ATL08 object (output of [ATL08_read()] function).
#' An S4 object of class [`ICESat2VegR::icesat2.atl08_dt-class`].
#' @param beam Character vector indicating beams to process
#' (e.g. "gt1l", "gt1r", "gt2l", "gt2r", "gt3l", "gt3r")
#'
#' @return Returns an S4 object of class [`ICESat2VegR::icesat2.atl03atl08_dt-class`]
#' containing the ATL08 computed photons attributes.
#'
#' @details These are the photons attributes extracted by default:
#' - `ph_segment_id`: Georeferenced segment id (20-m) associated with each photon.
#' - `lon_ph`: Longitude of each received photon. Computed from the ECEF Cartesian coordinates
#'             of the bounce point.
#' - `lat_ph`: Latitude of each received photon. Computed from the ECEF Cartesian coordinates
#'             of the bounce point.
#' - `h_ph`: Height of each received photon, relative to the WGS-84 ellipsoid including the
#'           geophysical corrections noted in section 6.0. Please note that neither the geoid,
#'           ocean tide nor the dynamic atmospheric corrections (DAC) are applied to the
#'           ellipsoidal heights.
#' - `quality_ph`: Indicates the quality of the associated photon. 0=nominal,
#'                 1=possible_afterpulse, 2=possible_impulse_response_effect, 3=possible_tep.
#'                 Use this flag in conjunction with `signal_conf_ph` to identify those photons
#'                 that are likely noise or likely signal.
#' - `solar_elevation`: Elevation of the sun above the horizon at the photon bounce point.
#' - `dist_ph_along`: Along-track distance of the photon from the beginning of the segment.
#' - `dist_ph_across`: Across-track distance of the photon from the center of the segment.
#' - `night_flag`: Flag indicating the data were acquired in night conditions: 0=day, 1=night.
#'                 Night flag is set when solar elevation is below 0.0 degrees.
#' - `classed_pc_indx`: Indices of photons tracking back to ATL03 that surface finding software
#'                      identified and used within the creation of the data products.
#' - `classed_pc_flag`: The L2B algorithm is run if this flag is set to 1 indicating data have
#'                      sufficient waveform fidelity for L2B to run.
#' - `ph_h`: Height of photon above interpolated ground surface.
#' - `d_flag`: Flag indicating whether DRAGANN labeled the photon as noise or signal.
#' - `delta_time`: Mid-segment GPS time in seconds past an epoch. The epoch is provided in
#'                 the metadata at the file level.
#' - `orbit_number`: Orbit number identifier to identify data from different orbits.
#' - `beam`: Beam identifier.
#' - `strong_beam`: Logical indicating if the beam is a strong beam.
#'
#' @seealso \url{https://icesat-2.gsfc.nasa.gov/sites/default/files/page_files/ICESat2_ATL08_ATBD_r006.pdf}
#'
#' @examples
#' # Specifying the path to ATL03 file
#' atl03_path <- system.file("extdata",
#'   "atl03_clip.h5",
#'   package = "ICESat2VegR"
#' )
#'
#' # Specifying the path to ATL08 file
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
#' # # Extracting ATL03 and ATL08 photons and heights
#' atl03_atl08_dt <- ATL03_ATL08_photons_attributes_dt_join(atl03_h5, atl08_h5)
#' head(atl03_atl08_dt)
#'
#' close(atl03_h5)
#' close(atl08_h5)
#' @export
ATL03_ATL08_photons_attributes_dt_join <- function(
    atl03_h5, atl08_h5,
    beam = c("gt1l", "gt1r", "gt2l", "gt2r", "gt3l", "gt3r")) {
  night_flag <-
    ph_segment_id <-
    classed_pc_indx <-
    . <-
    solar_elevation <- NA

  `:=` <- data.table::`:=`
  `.N` <- data.table::`.N`

  # Check file input
  if (!inherits(atl03_h5, "icesat2.atl03_h5")) {
    stop("atl03_h5 must be an object of class 'icesat2.atl03_h5' - output of [ATL03_read()] function ")
  }

  # Check file input
  if (!inherits(atl08_h5, "icesat2.atl08_h5")) {
    stop("atl08_h5 must be an object of class 'icesat2.atl08_h5' - output of [ATL08_read()] function ")
  }

  # Check beams to select
  beam_atl03 <- intersect(beam, atl03_h5$beams)
  beam_atl08 <- intersect(beam, atl08_h5$beams)

  beam <- intersect(beam_atl03, beam_atl08)

  photon.dt <- list()
  pb <- utils::txtProgressBar(min = 0, max = length(beam), style = 3)

  i_s <- 0
  for (i in beam) {
    atl03_beam_group <- atl03_h5[[i]]
    atl08_beam_group <- atl08_h5[[i]]
    i_s <- i_s + 0.25
    utils::setTxtProgressBar(pb, i_s)

    if (!atl03_beam_group$exists("geolocation/segment_ph_cnt")) {
      i_s <- i_s + 0.75
      utils::setTxtProgressBar(pb, i_s)
      next
    }

    segment_ph_cnt <- atl03_beam_group[["geolocation/segment_ph_cnt"]][]

    dataTableATL03Segs <- data.table::data.table(
      ph_segment_id = atl03_beam_group[["geolocation/segment_id"]][],
      ph_index_beg = atl03_beam_group[["geolocation/ph_index_beg"]][],
      segment_ph_cnt = segment_ph_cnt
    )

    dataTableATL08Photons <- data.table::data.table(
      ph_segment_id = atl08_beam_group[["signal_photons/ph_segment_id"]][],
      classed_pc_indx = atl08_beam_group[["signal_photons/classed_pc_indx"]][],
      classed_pc_flag = atl08_beam_group[["signal_photons/classed_pc_flag"]][],
      ph_h = atl08_beam_group[["signal_photons/ph_h"]][],
      d_flag = atl08_beam_group[["signal_photons/d_flag"]][],
      delta_time = atl08_beam_group[["signal_photons/delta_time"]][]
    )

    i_s <- i_s + 0.25
    utils::setTxtProgressBar(pb, i_s)

    n_segments <- atl03_beam_group[["geolocation/segment_length"]]$dims
    if (n_segments > 1) {
      segment_length <- c(0, cumsum(atl03_beam_group[["geolocation/segment_length"]][1:(n_segments - 1)]))
    } else {
      segment_length <- c(0)
    }
    segment_lengths <- rep(segment_length, segment_ph_cnt)

    segment_solar_elevation <- atl03_beam_group[["geolocation/solar_elevation"]][]
    ph_solar_elev <- rep(segment_solar_elevation, segment_ph_cnt)

    ph_segment_id <- rep(dataTableATL03Segs$ph_segment_id, dataTableATL03Segs$segment_ph_cnt)

    dataTableATL03Photons <- data.table::data.table(
      ph_segment_id = ph_segment_id,
      lon_ph = atl03_beam_group[["heights/lon_ph"]][],
      lat_ph = atl03_beam_group[["heights/lat_ph"]][],
      h_ph = atl03_beam_group[["heights/h_ph"]][],
      quality_ph = atl03_beam_group[["heights/quality_ph"]][],
      solar_elevation = ph_solar_elev,
      dist_ph_along = atl03_beam_group[["heights/dist_ph_along"]][] + segment_lengths,
      dist_ph_across = atl03_beam_group[["heights/dist_ph_across"]][]
    )

    i_s <- i_s + 0.25
    utils::setTxtProgressBar(pb, i_s)

    dataTableATL03Photons[, night_flag := as.integer(solar_elevation < 0)]
    unique_segment_id <- unique(dataTableATL08Photons$ph_segment_id)
    unique_segment_id <- intersect(unique_segment_id, dataTableATL03Segs$ph_segment_id)

    dataTableATL08Photons <- dataTableATL08Photons[ph_segment_id %in% unique_segment_id]
    dataTableATL03Photons <- dataTableATL03Photons[ph_segment_id %in% unique_segment_id]
    dataTableATL03Photons[, classed_pc_indx := seq_len(.N), by = ph_segment_id]
    
    data.table::setindex(dataTableATL03Photons, ph_segment_id, classed_pc_indx)
    data.table::setindex(dataTableATL08Photons, ph_segment_id, classed_pc_indx)
    ph_dt <- dataTableATL03Photons[dataTableATL08Photons, on = .(ph_segment_id, classed_pc_indx)]
    ph_dt$orbit_number <- atl03_h5[["orbit_info/orbit_number"]][]
    ph_dt[, beam := i]
    ph_dt$strong_beam <- i %in% atl03_h5$strong_beams


    photon.dt[[""]] <- ph_dt

    i_s <- i_s + 0.25
    utils::setTxtProgressBar(pb, i_s)
  }

  photon_dt <- data.table::rbindlist(photon.dt, fill = TRUE)

  prepend_class(photon_dt, "icesat2.atl03atl08_dt")
  close(pb)
  return(photon_dt)
}
