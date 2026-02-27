ATL03.photon.map <- list()

# Core geolocation
ATL03.photon.map[["delta_time"]]      <- "heights/delta_time"
ATL03.photon.map[["dist_ph_across"]]  <- "heights/dist_ph_across"
ATL03.photon.map[["dist_ph_along"]]  <- "heights/dist_ph_along"

# Photon ID
ATL03.photon.map[["pce_mframe_cnt"]]  <- "heights/pce_mframe_cnt"
ATL03.photon.map[["ph_id_channel"]]   <- "heights/ph_id_channel"
ATL03.photon.map[["ph_id_count"]]     <- "heights/ph_id_count"
ATL03.photon.map[["ph_id_pulse"]]     <- "heights/ph_id_pulse"

# Quality & signal
ATL03.photon.map[["quality_ph"]]      <- "heights/quality_ph"
ATL03.photon.map[["signal_class_ph"]] <- "heights/signal_class_ph"
ATL03.photon.map[["signal_conf_ph"]]  <- "heights/signal_conf_ph"
ATL03.photon.map[["weight_ph"]]       <- "heights/weight_ph"

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
#' @param attributes
#' Character vector of additional photon-level attributes to extract.
#' See details below for available options. Default is
#' `c("quality_ph", "dist_ph_along")`.
#'
#' @return
#' An S4 object of class [`data.table::data.table`] (with class
#' `icesat2.atl03_dt` prepended) containing photon-level ATL03
#' attributes.
#'
#' @details
#' The returned `data.table` always includes the following columns:
#'
#' - `beam`: Beam ID (e.g., `"gt2l"`).
#' - `strong_beam`: Logical indicating whether the beam is classified as
#'   strong for that orbit.
#' - `lon_ph`: Longitude of each received photon, computed from ECEF
#'   Cartesian coordinates of the bounce point (degrees_east).
#' - `lat_ph`: Latitude of each received photon, computed from ECEF
#'   Cartesian coordinates of the bounce point (degrees_north).
#' - `h_ph`: Height of each received photon, relative to the WGS-84
#'   ellipsoid, including the geophysical corrections described in the
#'   ATL03 ATBD. (Geoid, ocean tide, and dynamic atmosphere corrections
#'   are *not* applied.)
#' - `solar_elevation`: Solar elevation interpolated from segment-level
#'   `geolocation/solar_elevation` to photon level.
#'
#' Additional photon-level attributes can be requested using the
#' `attributes` argument. These must correspond to names defined in
#' `ATL03.photon.map`. Available optional attributes currently include:
#'
#' - `delta_time`: Photon transmit time in seconds since
#'   2018-01-01 (ATLAS SDP epoch).
#' - `dist_ph_across`: Across-track distance of the photon projected
#'   to the reference ellipsoid (meters).
#' - `pce_mframe_cnt`: Major frame counter (part of photon ID).
#' - `ph_id_channel`: Channel number assigned to the photon event.
#' - `ph_id_count`: Photon event counter (part of photon ID).
#' - `ph_id_pulse`: Laser pulse counter (part of photon ID).
#' - `quality_ph`: Photon quality flag indicating nominal, saturation,
#'   or noise conditions.
#' - `signal_class_ph`: Experimental photon classification flag based
#'   on photon weights.
#' - `signal_conf_ph`: Signal confidence level per surface type
#'   (5 x N array in original product).
#' - `weight_ph`: Relative photon weight (0-65535), proportional to
#'   photon density.
#'
#' If a requested attribute is not present in a given beam, the
#' corresponding column is filled with `NA`.
#'
#' @seealso
#' \url{https://nsidc.org/sites/default/files/documents/technical-reference/icesat2_atl03_data_dict_v007.pdf}
#'
#' @examples
#' atl03_path <- system.file(
#'   "extdata", "atl03_clip.h5",
#'   package = "ICESat2VegR"
#' )
#'
#' atl03_h5 <- ATL03_read(atl03_path)
#'
#' atl03_photons_dt <- ATL03_photons_attributes_dt(
#'  atl03_h5,
#'  beam = c("gt1r"),
#'  attributes = c("quality_ph", "dist_ph_along", "ph_id_count")
#' )
#' head(atl03_photons_dt)
#'
#' close(atl03_h5)
#'
#'
#' @export
ATL03_photons_attributes_dt <- function(
  atl03_h5,
  beam = c("gt1l", "gt1r", "gt2l", "gt2r", "gt3l", "gt3r"),
  attributes = c(
    "quality_ph",
    "dist_ph_along",
  )
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
    if (n_segments > 1 && "dist_ph_along" %in% attributes) {
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
      segment_solar_elevation <- atl03_h5[[
          paste0(i, "/geolocation/solar_elevation")
        ]][]
    } else {
      segment_solar_elevation <- rep(NA_real_, n_segments)
    }
    ph_solar_elev <- rep(segment_solar_elevation, segment_ph_cnt)

    # --------------------------------------------------------------
    # Photon-level quality
    # --------------------------------------------------------------
    heights_group <- paste0(i, "/heights/")
    datasets_heights <- atl03_h5[[heights_group]]$dt_datasets()$name
    
    # number of photons (reference length)
    n_ph <- length(atl03_h5[[paste0(i, "/heights/lon_ph")]][])

    # container for extracted attributes
    extra_cols <- list()

    # Validate requested attributes
    attributes <- intersect(attributes, names(ATL03.photon.map))

    for (attr in attributes) {

      dataset_name <- basename(ATL03.photon.map[[attr]])

      if (dataset_name %in% datasets_heights) {

        values <- atl03_h5[[paste0(i, "/", ATL03.photon.map[[attr]])]][]

        # Special handling for dist_ph_along
        if (attr == "dist_ph_along") {
          values <- values + segment_lengths
        }

      } else {
        # missing dataset -> NA vector
        values <- rep(NA, n_ph)
      }

      extra_cols[[attr]] <- values
    }



    # --------------------------------------------------------------
    # Assemble photon-level data.table for this beam
    # --------------------------------------------------------------
    dt_atl03_photons <- data.table::data.table(
      beam            = i,
      strong_beam     = i %in% atl03_h5$strong_beams,
      lon_ph          = atl03_h5[[paste0(i, "/heights/lon_ph")]][],
      lat_ph          = atl03_h5[[paste0(i, "/heights/lat_ph")]][],
      h_ph            = atl03_h5[[paste0(i, "/heights/h_ph")]][],
      solar_elevation = ph_solar_elev
    )

    # Append optional attributes
    for (nm in names(extra_cols)) {
      dt_atl03_photons[[nm]] <- extra_cols[[nm]]
    }

    photon_dt <- data.table::rbindlist(
      list(photon_dt, dt_atl03_photons),
      fill = TRUE
    )

    utils::setTxtProgressBar(pb, i_s)
  }

  close(pb)

  prepend_class(photon_dt, "icesat2.atl03_dt")
  photon_dt
}
