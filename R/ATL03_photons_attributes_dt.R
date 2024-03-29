#' ATL03 photons attributes
#'
#' @description This function extracts photons attributes from ICESat-2 ATL03 data
#'
#' @usage ATL03_photons_attributes_dt(ATL03_h5, beam)
#'
#' @param ATL03_h5 A ICESat-2 ATL03 object (output of [ATL03_read()] function).
#' An S4 object of class [`ICESat2VegR::icesat2.atl03_dt-class`].
#' @param beam Character vector indicating beams to process (e.g. "gt1l", "gt1r", "gt2l", "gt2r", "gt3l", "gt3r")
#'
#' @return Returns an S4 object of class [`data.table::data.table-class`]
#' containing the ATL03 photons attributes.
#'
#' @details These are the photons attributes extracted:
#' - \emph{lon_ph} Longitude of each received photon.
#' Computed from the ECEF Cartesian coordinates of the bounce point.
#' - \emph{lat_ph} Latitude of each received photon.
#' Computed from the ECEF Cartesian coordinates of the bounce point.
#' - \emph{lat_ph} Latitude of each received photon.
#' Computed from the ECEF Cartesian coordinates of the bounce point.
#' Height of each received photon, relative to the WGS-84 ellipsoid including
#' the geophysical corrections noted in section 6.0. Please note that
#' neither the geoid, ocean tide nor the dynamic atmospheric corrections (DAC) are applied to the ellipsoidal heights.
#' - \emph{quality_ph} Indicates the quality of the associated photon.
#' 0=nominal, 1=possible_afterpulse, 2=possible_impulse_response_
#' effect, 3=possible_tep. Use this flag in conjunction with signal_conf_ph
#' to identify those photons that are likely noise or likely signal
#' - \emph{night_flag}  Flag indicating the data were acquired in night
#' conditions: 0=day, 1=night. Night flag is set when solar elevation is below 0.0 degrees.
#'
#'
#' @seealso \url{https://icesat-2.gsfc.nasa.gov/sites/default/files/page_files/ICESat2_ATL03_ATBD_r006.pdf}
#'
#' @examples
#'
#' # Specifying the path to ATL03 file (zip file)
#' atl03_path <- system.file("extdata",
#'   "atl03_clip.h5",
#'   package = "ICESat2VegR"
#' )
#'
#' # Reading ATL03 data (h5 file)
#' atl03_h5 <- ATL03_read(atl03_path = atl03_path)
#'
#' # Extracting ATL03 photons attributes
#' atl03_photons_dt <- ATL03_photons_attributes_dt(atl03_h5 = atl03_h5)
#'
#' head(atl03_photons_dt)
#' close(atl03_h5)
#' @export
ATL03_photons_attributes_dt <- function(atl03_h5,
                                        beam = c("gt1l", "gt1r", "gt2l", "gt2r", "gt3l", "gt3r")) {
  # Check file input
  if (!inherits(atl03_h5, "icesat2.atl03_h5")) {
    stop("atl03_h5 must be an object of class 'icesat2.atl03_h5' - output of [ATL03_read()] function ")
  }

  # Check beams to select
  beam <- intersect(beam, atl03_h5$beams)

  photon.dt <- data.table::data.table()

  pb <- utils::txtProgressBar(min = 0, max = length(beam), style = 3)

  i_s <- 0

  for (i in beam) {
    i_s <- i_s + 1
    message(i)

    if (!atl03_h5$exists(paste0(i, "/geolocation/segment_length"))) {
      next
    }
    n_segments <- atl03_h5[[paste0(i, "/geolocation/segment_length")]]$dims
    if (n_segments == 0) {
      next
    }
    segment_ph_cnt <- atl03_h5[[paste0(i, "/geolocation/segment_ph_cnt")]][]
    segment_length <- c(0, cumsum(atl03_h5[[paste0(i, "/geolocation/segment_length")]][1:(n_segments - 1)]))
    segment_lengths <- rep(segment_length, segment_ph_cnt)

    if ("solar_elevation" %in% atl03_h5[[paste0(i, "/geolocation/")]]$dt_datasets()$name) {
      segment_solar_elevation <- atl03_h5[[paste0(i, "/geolocation/solar_elevation")]][]
    } else {
      segment_solar_elevation <- rep(NA, n_segments)
    }
    ph_solar_elev <- rep(segment_solar_elevation, segment_ph_cnt)


    if ("quality_ph" %in% atl03_h5[[paste0(i, "/heights/")]]$dt_datasets()$name) {
      quality_ph <- atl03_h5[[paste0(i, "/heights/quality_ph")]][]
    } else {
      quality_ph <- rep(NA, length(atl03_h5[[paste0(i, "/heights/lon_ph")]][]))
    }


    if ("dist_ph_along" %in% atl03_h5[[paste0(i, "/heights/")]]$dt_datasets()$name) {
      dist_ph_along <- atl03_h5[[paste0(i, "/heights/dist_ph_along")]][] + segment_lengths
    } else {
      dist_ph_along <- rep(NA, length(atl03_h5[[paste0(i, "/heights/lon_ph")]][] + segment_lengths))
    }


    dataTableATL03Photons <- data.table::data.table(cbind(
      lon_ph = atl03_h5[[paste0(i, "/heights/lon_ph")]][],
      lat_ph = atl03_h5[[paste0(i, "/heights/lat_ph")]][],
      h_ph = atl03_h5[[paste0(i, "/heights/h_ph")]][],
      quality_ph = quality_ph,
      solar_elevation = ph_solar_elev,
      dist_ph_along = dist_ph_along
    ))

    photon_dt <- data.table::rbindlist(list(photon.dt, dataTableATL03Photons), fill = TRUE)
    utils::setTxtProgressBar(pb, i_s)
  }

  prepend_class(photon_dt, "icesat2.atl03_dt")


  close(pb)

  return(photon_dt)
}
