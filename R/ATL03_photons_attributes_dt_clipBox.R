#' Clip ATL03 photons by Coordinates
#'
#' @description This function clips ATL03 photons attributes within a given bounding coordinates
#'
#' @usage ATL03_photons_attributes_dt_clipBox(atl03_photons_dt, xmin, xmax, ymin, ymax)
#'
#' @param atl03_photons_dt A atl03_photons_dt object (output of [atl03_photons_attributes_dt()] function).
#' An S4 object of class [ICESat2VegR::icesat2.atl03_dt]
#' @param lower_left_lon Numeric. West longitude (x) coordinate of bounding rectangle, in decimal degrees.
#' @param upper_right_lon Numeric. East longitude (x) coordinate of bounding rectangle, in decimal degrees.
#' @param lower_left_lat Numeric. South latitude (y) coordinate of bounding rectangle, in decimal degrees.
#' @param upper_right_lat Numeric. North latitude (y) coordinate of bounding rectangle, in decimal degrees.
#'
#' @return Returns an S4 object of class [ICESat2VegR::icesat2.atl03_dt]
#' containing the ATL03 photons attributes.
#'
#' @seealso \url{https://icesat-2.gsfc.nasa.gov/sites/default/files/page_files/ICESat2_atl03_ATBD_r006.pdf}
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
#' # Extracting ATL03 photons attributes
#' atl03_photons_dt <- ATL03_photons_attributes_dt(atl03_h5 = atl03_h5)
#'
#' # Bounding rectangle coordinates
#' lower_left_lon <- -103.7604
#' lower_left_lat <- 59.4672
#' upper_right_lon <- -103.7600
#' upper_right_lat <- 59.4680
#'
#' # Clipping ATL08-derived canopy metrics by boundary box extent
#' atl03_photons_dt_clip <- ATL03_photons_attributes_dt_clipBox(
#'   atl03_photons_dt,
#'   lower_left_lon,
#'   upper_right_lon,
#'   upper_right_lat,
#'   lower_left_lat
#' )
#'
#' head(atl03_photons_dt_clip)
#'
#' close(atl03_h5)
#' @import hdf5r stats
#' @export
ATL03_photons_attributes_dt_clipBox <- function(atl03_photons_dt,
                                                lower_left_lon,
                                                upper_right_lon,
                                                upper_right_lat,
                                                lower_left_lat) {
  if (!inherits(atl03_photons_dt, "icesat2.atl03_dt")) {
    stop("atl03_photons_dt needs to be an object of class 'icesat2.atl03_dt' ")
  }

  if (any(is.na(atl03_photons_dt))) {
    atl03_photons_dt <- na.omit(atl03_photons_dt)
  }


  # xmin ymin xmax ymax
  mask <-
    atl03_photons_dt$lon_ph >= lower_left_lon &
      atl03_photons_dt$lon_ph <= upper_right_lon &
      atl03_photons_dt$lat_ph >= lower_left_lat &
      atl03_photons_dt$lat_ph <= upper_right_lat

  mask[!stats::complete.cases(mask)] <- FALSE
  mask <- (seq_along(atl03_photons_dt$lat_ph))[mask]
  newFile <- atl03_photons_dt[mask, ]

  # newFile<- new("gedi.level1b.dt", dt = level1bdt[mask,])
  return(newFile)
}
