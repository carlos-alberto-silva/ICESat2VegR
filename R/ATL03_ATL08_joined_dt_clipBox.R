#' Clip Joined ATL03 and ATL08 by Bounding Box
#'
#' @description This function clips joined ATL03 and ATL08 photon attributes table within a given Bounding Box
#'
#' @usage ATL03_ATL08_joined_dt_clipBox(atl03_atl08_dt, xmin, xmax, ymin, ymax)
#'
#' @param atl03_atl08_dt  An S4 object of class [rICESat2Veg::icesat2.atl03atl08_dt] containing ATL03 and ATL08 data
#' (output of [rICESat2Veg::ATL03_ATL08join()] function).
#' @param xmin Numeric. West longitude (x) coordinate of bounding rectangle, in decimal degrees.
#' @param xmax Numeric. East longitude (x) coordinate of bounding rectangle, in decimal degrees.
#' @param ymin Numeric. South latitude (y) coordinate of bounding rectangle, in decimal degrees.
#' @param ymax Numeric. North latitude (y) coordinate of bounding rectangle, in decimal degrees.
#'
#' @return Returns an S4 object of class [rICESat2Veg::icesat2.atl08_dt]
#' containing the ATL08 Canopy Height Metrics.
#'
#'@seealso \url{https://icesat-2.gsfc.nasa.gov/sites/default/files/page_files/ICESat2_ATL08_ATBD_r006.pdf}
#'
#' @examples
#'# Specifying the path to ATL03 and ATL08 file (zip file)
#'outdir = tempdir()
#'atl03_zip <- system.file("extdata",
#'                   "ATL03_20220401221822_01501506_005_01.zip",
#'                   package="rICESat2Veg")
#'
#'atl08_zip <- system.file("extdata",
#'                   "ATL08_20220401221822_01501506_005_01.zip",
#'                   package="rICESat2Veg")
#'
#'# Unzipping ATL03 file
#'atl03_path <- unzip(atl03_zip,exdir = outdir)
#'
#'# Unzipping ATL08 file
#'atl08_path <- unzip(atl08_zip,exdir = outdir)
#'
#'# Reading ATL03 data (h5 file)
#atl03_h5<-ATL08read(atl03_path=atl03_path)
#'
#'# Reading ATL08 data (h5 file)
#atl08_h5<-ATL08read(atl08_path=atl08_path)
#'
#'# # Extracting ATL03 and ATL08 photons and heights
#'atl03_08_dt<-ATL03_ATL08join(atl03_h5,atl08_h5)
#'head(atl03_08_dt)
#'
#' # Bounding rectangle coordinates
#' xmin <- -107.7
#' xmax <- -106.5
#' ymin <- 32.75
#' ymax <- 42.75
#'
#' # Clipping ATL08-derived canopy metrics by boundary box extent
#'atl03_08_dt_clip <- ATL03_ATL08_joined_dt_clipBox(atl03_08_dt, xmin, xmax, ymin, ymax)
#'
#'close(atl03_h5)
#'close(atl08_h5)
#'@import hdf5r stats
#'@export
ATL03_ATL08_joined_dt_clipBox <- function(atl03_atl08_dt, xmin, xmax, ymin, ymax) {

  if (!class(atl03_atl08_dt)[1]=="icesat2.atl03atl08_dt"){
    stop("atl03_atl08_dt needs to be an object of class 'icesat2.at03atl08_dt' ")
  }

  if (any(is.na(atl03_atl08_dt@dt))) {
    atl03_atl08_dt<-na.omit(atl03_atl08_dt@dt)
  } else {

    atl03_atl08_dt<-atl03_atl08_dt@dt
  }


  # xmin ymin xmax ymax
  mask <-
    atl03_atl08_dt$lon_ph >= xmin &
    atl03_atl08_dt$lon_ph <= xmax &
    atl03_atl08_dt$lat_ph >= ymin &
    atl03_atl08_dt$lat_ph <= ymax

  mask[!stats::complete.cases(mask)] <- FALSE
  mask <- (seq_along(atl03_atl08_dt$lat_ph))[mask]
  newFile <- atl03_atl08_dt[mask, ]

  newFile<- new("icesat2.atl08_dt", dt = newFile)

  # newFile<- new("gedi.level1b.dt", dt = level1bdt[mask,])
  return(newFile)
}
