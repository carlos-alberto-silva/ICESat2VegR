#' ATL03 photon attributes
#'
#' @description
#' Extract photon-level attributes from ICESat-2 ATL03 data.
#'
#' @param atl03_h5
#'   An ICESat-2 ATL03 object (output of [ATL03_read()]), i.e.
#'   an S4 object of class [`ICESat2VegR::icesat2.atl03_h5`].
#' @param beam
#'   Character vector indicating beams to process
#'   (e.g. "gt1l", "gt1r", "gt2l", "gt2r", "gt3l", "gt3r").
#'
#' @return
#' An S4 object of class [`data.table::data.table`] (with class
#' `"icesat2.atl03_dt"` prepended) containing photon-level ATL03
#' attributes.
#'
#' @details
#' Attributes extracted include:
#'
#' - `lon_ph`: Longitude of each received photon, computed from ECEF
#'   Cartesian coordinates of the bounce point.
#' - `lat_ph`: Latitude of each received photon, computed from ECEF
#'   Cartesian coordinates of the bounce point.
#' - `h_ph`: Height of each received photon, relative to the WGS-84
#'   ellipsoid, including the geophysical corrections noted in the
#'   ATL03 ATBD. (Geoid, ocean tide and DAC are *not* applied.)
#' - `quality_ph`: Photon quality flag
#'   (0 = nominal, 1 = possible_afterpulse, 2 = possible_impulse_response_effect,
#'   3 = possible_tep). Use together with `signal_conf_ph` to identify
#'   likely noise vs likely signal.
#' - `solar_elevation`: Solar elevation interpolated from segment-level
#'   `geolocation/solar_elevation` to each photon.
#' - `dist_ph_along`: Along-track distance for each photon (segment
#'   cumulative length + `heights/dist_ph_along` when available).
#' - `beam`: Beam ID (e.g. "gt2l").
#' - `strong_beam`: Logical indicating whether the beam is classified as
#'   strong for that orbit.
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
#' atl03_photons_dt <- ATL03_photons_attributes_dt(atl03_h5)
#' head(atl03_photons_dt)
#'
#' close(atl03_h5)
#' }
#'
#' @export
ATL03_photons_attributes_dt <- function(
  atl03_h5,
  beam = c("gt1l", "gt1r", "gt2l", "gt2r", "gt3l", "gt3r")
) {
  # ----------------------
  # Basic input checks
  # ----------------------
  if (!inherits(atl03_h5, "icesat2.atl03_h5")) {
    stop("atl03_h5 must be an object of class 'icesat2.atl03_h5' (output of ATL03_read()).")
  }

  # Restrict to beams actually present in the file
  beam <- intersect(beam, atl03_h5$beams)
  if (length(beam) == 0) {
    stop("None of the requested beams are present in 'atl03_h5'.")
  }

  photon_dt <- data.table::data.table()

  pb <- utils::txtProgressBar(
    min   = 0,
    max   = length(beam),
    style = 3,
    file  = stderr()
  )

  i_s <- 0

  for (i in beam) {
    i_s <- i_s + 1

    # --------------------------------------------------------------
    # Check that we have segment info; otherwise skip the beam
    # --------------------------------------------------------------
    seg_len_name <- paste0(i, "/geolocation/segment_length")
    if (!atl03_h5$exists(seg_len_name)) {
      utils::setTxtProgressBar(pb, i_s)
      next
    }

    n_segments <- atl03_h5[[seg_len_name]]$dims
    if (n_segments == 0) {
      utils::setTxtProgressBar(pb, i_s)
      next
    }

    seg_ph_cnt_name <- paste0(i, "/geolocation/segment_ph_cnt")
    segment_ph_cnt  <- atl03_h5[[seg_ph_cnt_name]][]

    if (length(segment_ph_cnt) == 0) {
      utils::setTxtProgressBar(pb, i_s)
      next
    }

    # --------------------------------------------------------------
    # Segment-length expansion for dist_ph_along
    # --------------------------------------------------------------
    if (n_segments > 1) {
      seg_len <- atl03_h5[[seg_len_name]][1:(n_segments - 1)]
      segment_length_offset <- c(0, cumsum(seg_len))
      segment_lengths       <- rep(segment_length_offset, segment_ph_cnt)
    } else {
      segment_lengths <- 0
    }

    # --------------------------------------------------------------
    # Solar elevation (segment -> photon)
    # --------------------------------------------------------------
    geol_group <- paste0(i, "/geolocation/")
    datasets_geol <- atl03_h5[[geol_group]]$dt_datasets()$name

    if ("solar_elevation" %in% datasets_geol) {
      segment_solar_elevation <- atl03_h5[[paste0(i, "/geolocation/solar_elevation")]][]
    } else {
      segment_solar_elevation <- rep(NA_real_, n_segments)
    }
    ph_solar_elev <- rep(segment_solar_elevation, segment_ph_cnt)

    # --------------------------------------------------------------
    # Photon-level quality
    # --------------------------------------------------------------
    heights_group <- paste0(i, "/heights/")
    datasets_heights <- atl03_h5[[heights_group]]$dt_datasets()$name

    if ("quality_ph" %in% datasets_heights) {
      quality_ph <- atl03_h5[[paste0(i, "/heights/quality_ph")]][]
    } else {
      n_ph <- length(atl03_h5[[paste0(i, "/heights/lon_ph")]][])
      quality_ph <- rep(NA_integer_, n_ph)
    }

    # --------------------------------------------------------------
    # dist_ph_along
    # --------------------------------------------------------------
    if ("dist_ph_along" %in% datasets_heights) {
      dist_ph_along <- atl03_h5[[paste0(i, "/heights/dist_ph_along")]][] +
        segment_lengths
    } else {
      n_ph <- length(atl03_h5[[paste0(i, "/heights/lon_ph")]][])
      dist_ph_along <- rep(NA_real_, n_ph)
    }

    # --------------------------------------------------------------
    # Assemble photon-level data.table for this beam
    # --------------------------------------------------------------
    dataTableATL03Photons <- data.table::data.table(
      beam           = i,
      strong_beam    = i %in% atl03_h5$strong_beams,
      lon_ph         = atl03_h5[[paste0(i, "/heights/lon_ph")]][],
      lat_ph         = atl03_h5[[paste0(i, "/heights/lat_ph")]][],
      h_ph           = atl03_h5[[paste0(i, "/heights/h_ph")]][],
      quality_ph     = quality_ph,
      solar_elevation = ph_solar_elev,
      dist_ph_along  = dist_ph_along
    )

    photon_dt <- data.table::rbindlist(
      list(photon_dt, dataTableATL03Photons),
      fill = TRUE
    )

    utils::setTxtProgressBar(pb, i_s)
  }

  close(pb)

  prepend_class(photon_dt, c("icesat2.atl03_dt","data.table","data.frame"))
  photon_dt
}
