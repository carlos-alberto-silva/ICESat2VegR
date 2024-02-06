#' Extract reference ground track from ATL03 segments
#'
#' @description This function extracts reference ground track from ICESat-2 ATL03 data
#' and writes it as a GDAL vector format.
#'
#' @param ATL03_h5 A ICESat-2 ATL03 object (output of [ATL03_read()] function).
#' An S4 object of class [ICESat2VegR::icesat2.atl03_dt].
#' @param output Character vector. The output vector file. The GDAL vector format
#' will be inferred by the file extension using [`terra::writeVector()`]
#' @param overwrite logical input to control if the output vector file should be overwritten.
#' Default FALSE.
#'
#' @return Returns the lines vector representing the ground track
#'
#' @details This function will use the reference photons from the segments as reference
#' for deriving the ground tracks. The begining and end of the lines are interpolated from
#' the information regarding the position of the reference photon within the segment and the
#' segment length.
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
#' atl03_photons_dt <- ATL03_rgt_extract(atl03_h5 = atl03_h5)
#' head(ATL03_photons_dt)
#'
#' close(ATL03_h5)
#' @include class.icesat2.R
#' @importFrom data.table data.table rbindlist
#' @importFrom terra writeVector
#' @export
ATL03_rgt_extract <- function(atl03_h5, output, beam = c("gt1r", "gt2r", "gt1l", "gt3l", "gt2l", "gt3r"), overwrite = FALSE) {
  beams <- beam

  all_pts <- list()
  all_attr <- list()

  pb <- utils::txtProgressBar(min = 0, max = length(beams), style = 3)
  i_s <- 0

  for (beam in beams) {
    i_s <- i_s + 1

    ph_index_beg <- atl03_h5[[paste0(beam, "/geolocation/ph_index_beg")]][]
    photon_mask <- ph_index_beg > 0

    n_mask <- sum(photon_mask)

    reference_photon_lat <- atl03_h5[[paste0(beam, "/geolocation/reference_photon_lat")]][photon_mask]
    reference_photon_lon <- atl03_h5[[paste0(beam, "/geolocation/reference_photon_lon")]][photon_mask]
    segment_id <- atl03_h5[[paste0(beam, "/geolocation/segment_id")]][photon_mask]

    reference_photon_index <- atl03_h5[[paste0(beam, "/geolocation/reference_photon_index")]][]
    dist_ph_along <- atl03_h5[[
      paste0(beam, "/heights/dist_ph_along")
    ]][
      (ph_index_beg + reference_photon_index - 1)[photon_mask]
    ]

    segment_length <- atl03_h5[[paste0(beam, "/geolocation/segment_length")]][]
    segment_cumsum <- cumsum(segment_length)[photon_mask]

    distance_offset <- (diff(segment_cumsum) - segment_length[photon_mask][-n_mask])
    dist_head <- dist_ph_along
    dist_tail <- segment_length[photon_mask] - dist_ph_along
    distTotal <- distance_offset + dist_tail[-n_mask] + dist_head[-1]


    center <- data.table::data.table(
      lat = reference_photon_lat,
      lon = reference_photon_lon
    )

    center_diff <- center[, lapply(.SD, diff)]


    tail <- data.table::data.table(
      lat = numeric(length = n_mask),
      lon = numeric(length = n_mask)
    )

    head <- data.table::data.table(
      lat = numeric(length = n_mask),
      lon = numeric(length = n_mask)
    )

    #
    # The head and tail coordinates are computed as an offset from the
    # adjacent segment reference photons, adjusted by the distance
    # in meters from ATL03 data
    #
    tail[-n_mask] <- center[-n_mask] + (center_diff * dist_tail[-n_mask] / distTotal)
    head[-1] <- center[-1] - (center_diff * dist_head[-1] / distTotal)
    head[1] <- center[1] - dist_head[1] * ((tail[1] - center[1]) / dist_tail[1])

    tail[n_mask] <- center[n_mask] + ((center[n_mask] - head[n_mask]) * dist_tail[n_mask] / dist_head[n_mask])

    pts <- data.table::rbindlist(list(head[, list(object = .I, .SD)], tail[, list(object = .I, .SD)]))
    pts <- pts[order(object), list(
      object,
      x = .SD.lon,
      y = .SD.lat
    )]

    all_pts[[beam]] <- pts
    all_attr[[beam]] <- data.table::data.table(
      segment_id = segment_id,
      reference_photon_index = reference_photon_index[photon_mask],
      segment_length = segment_length[photon_mask],
      beam = beam
    )

    utils::setTxtProgressBar(pb, i_s)
  }

  all_pts2 <- data.table::rbindlist(all_pts)
  all_attr2 <- data.table::rbindlist(all_attr)

  v <- terra::vect(
    as.matrix(all_pts2),
    "lines",
    crs = "epsg:4326",
    atts = all_attr2
  )


  terra::writeVector(v, output, overwrite = TRUE)

  return(v)
}
