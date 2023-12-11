#' Clip ATL08 Canopy Height Metrics by Coordinates
#'
#' @description This function clips ATL08 Canopy Height Metrics within a given bounding coordinates
#'
#' @usage ATL08_canopy_dt_clipBox(atl08_canopy_dt, xmin, xmax, ymin, ymax)
#'
#' @param atl08_canopy_dt A atl08_canopy_dt object (output of [ATL08_canopy_attributes()] function).
#' An S4 object of class [rICESat2Veg::icesat2.atl08_dt]
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
#'# Specifying the path to ATL08 file (zip file)
#'outdir = tempdir()
#'atl08_zip <- system.file("extdata",
#'                   "ATL08_20220401221822_01501506_005_01.zip",
#'                   package="rICESat2Veg")
#'
#'# Unzipping ATL08 file
#'atl08_path <- unzip(atl08_zip,exdir = outdir)
#'
#'# Reading ATL08 data (h5 file)
#atl08_h5<-readATL08(ATL08path=atl08_path)
#'
#'# Extracting ATL08-derived Canopy Metrics
#'atl08_canopy_dt<-ATL08_canopy_attributes(atl08_h5=atl08_h5)
#'
#' # Bounding rectangle coordinates
#' xmin <- -107.7
#' xmax <- -106.5
#' ymin <- 32.75
#' ymax <- 42.75
#'
#' # Clipping ATL08-derived canopy metrics by boundary box extent
#'atl08_canopy_dt_clip <- ATL08_canopy_dt_clipBox(atl08_canopy_dt, xmin, xmax, ymin, ymax)
#'
#'close(level2a)
#'@import hdf5r stats
#'@export
ATL08_canopy_dt_clipBox <- function(atl08_canopy_dt, xmin, xmax, ymin, ymax) {

  if (!class(atl08_canopy_dt)[1]=="icesat2.atl08_dt"){
    stop("atl08_canopy_dt needs to be an object of class 'icesat2.atl08_dt' ")
  }

  if (any(is.na(atl08_canopy_dt@dt))) {
    atl08_canopy_dt<-na.omit(atl08_canopy_dt@dt)
  } else {

    atl08_canopy_dt<-atl08_canopy_dt@dt
  }


  # xmin ymin xmax ymax
  mask <-
    atl08_canopy_dt$longitude >= xmin &
    atl08_canopy_dt$longitude <= xmax &
    atl08_canopy_dt$latitude >= ymin &
    atl08_canopy_dt$latitude <= ymax

  mask[!stats::complete.cases(mask)] <- FALSE
  mask <- (seq_along(atl08_canopy_dt$latitude))[mask]
  newFile <- atl08_canopy_dt[mask, ]

  newFile<- new("icesat2.atl08_dt", dt = newFile)

  # newFile<- new("gedi.level1b.dt", dt = level1bdt[mask,])
  return(newFile)
}
