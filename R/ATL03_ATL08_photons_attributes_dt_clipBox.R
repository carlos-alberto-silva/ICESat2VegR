#' Clip joined ATL03 and ATL08 photons by Bounding Box
#'
#' @description This function clips joined ATL03 and ATL08 photon attributes within a given Bounding Box
#'
#' @usage ATL03_ATL08_photons_attributes_dt_clipBox(atl03_atl08_dt, xmin, xmax, ymin, ymax)
#'
#' @param atl03_atl08_dt  An S4 object of class [ICESat2VegR::icesat2.atl03atl08_dt] containing ATL03 and ATL08 data
#' (output of [ICESat2VegR::ATL03_ATL08_photons_attributes_dt_join()] function).
#' @param lower_left_lon Numeric. West longitude (x) coordinate of bounding rectangle, in decimal degrees.
#' @param upper_right_lon Numeric. East longitude (x) coordinate of bounding rectangle, in decimal degrees.
#' @param lower_left_lat Numeric. South latitude (y) coordinate of bounding rectangle, in decimal degrees.
#' @param upper_right_lat Numeric. North latitude (y) coordinate of bounding rectangle, in decimal degrees.
#'
#' @return Returns an S4 object of class [ICESat2VegR::icesat2.atl03atl08_dt]
#' containing a subset of the ATL03 and ATL08 photon attributes.
#'
#'@seealso \url{https://icesat-2.gsfc.nasa.gov/sites/default/files/page_files/ICESat2_ATL03_ATBD_r006.pdf}
#'@seealso \url{https://icesat-2.gsfc.nasa.gov/sites/default/files/page_files/ICESat2_ATL08_ATBD_r006.pdf}
#'
#' @examples
#'# Specifying the path to ATL03 and ATL08 file (zip file)
#'outdir = tempdir()
#'atl03_zip <- system.file("extdata",
#'                   "ATL03_20220401221822_01501506_005_01.zip",
#'                   package="ICESat2VegR")
#'
#'atl08_zip <- system.file("extdata",
#'                   "ATL08_20220401221822_01501506_005_01.zip",
#'                   package="ICESat2VegR")
#'
#'# Unzipping ATL03 file
#'atl03_path <- unzip(atl03_zip,exdir = outdir)
#'
#'# Unzipping ATL08 file
#'atl08_path <- unzip(atl08_zip,exdir = outdir)
#'
#'# Reading ATL03 data (h5 file)
#atl03_h5<-ATL08_read(atl03_path=atl03_path)
#'
#'# Reading ATL08 data (h5 file)
#atl08_h5<-ATL08_read(atl08_path=atl08_path)
#'
#'# Joining ATL03 and ATL08 photons and heights
#'atl03_atl08_dt<-ATL03_ATL08_photons_attributes_dt_join(atl03_h5,atl08_h5)
#'head(atl03_atl08_dt)
#'
#'# Bounding rectangle coordinates
#'lower_left_lon=-107.7
#'lower_left_lat=42.75
#'upper_right_lon=-106.5
#'upper_right_lat=32.75
#'
#'# Clipping ATL08-derived canopy metrics by boundary box extent
#'atl03_atl08_dt_clip <- ATL03_ATL08_photons_attributes_dt_clipBox(atl03_atl08_dt,
#'                                                                lower_left_lon,
#'                                                                upper_right_lon,
#'                                                                upper_right_lat,
#'                                                                lower_left_lat)
#'head(atl03_atl08_dt_clip)
#'
#'close(atl03_h5)
#'close(atl08_h5)
#'@import hdf5r stats
#'@export
ATL03_ATL08_photons_attributes_dt_clipBox <- function(atl03_atl08_dt,
                                                      lower_left_lon,
                                                      upper_right_lon,
                                                      upper_right_lat,
                                                      lower_left_lat) {

  if (!inherits(atl03_atl08_dt, "icesat2.atl03atl08_dt")){
    stop("atl03_atl08_dt needs to be an object of class 'icesat2.at03atl08_dt' ")
  }

  if (any(is.na(atl03_atl08_dt))) {
    atl03_atl08_dt<-na.omit(atl03_atl08_dt)
  }

  # xmin ymin xmax ymax
  mask <-
    atl03_atl08_dt$lon_ph >= lower_left_lon &
    atl03_atl08_dt$lon_ph <= upper_right_lon &
    atl03_atl08_dt$lat_ph >= upper_right_lat &
    atl03_atl08_dt$lat_ph <= lower_left_lat

  mask[!stats::complete.cases(mask)] <- FALSE
  mask <- (seq_along(atl03_atl08_dt$lat_ph))[mask]
  newFile <- atl03_atl08_dt[mask, ]

  # newFile<- new("gedi.level1b.dt", dt = level1bdt[mask,])
  return(newFile)
}
