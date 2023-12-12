#' Clip ATL08 Terrain attributes by Coordinates
#'
#' @description This function clips ATL08 terrain attributes within a given bounding coordinates
#'
#' @usage ATL08_terrain_dt_clipBox(atl08_terrain_dt, xmin, xmax, ymin, ymax)
#'
#' @param atl08_terrain_dt A atl08_terrain_dt object (output of [rICESat2Veg::ATL08_terrain_attributes_dt()] function).
#' An S4 object of class [rICESat2Veg::icesat2.atl08_dt]
#' @param xmin Numeric. West longitude (x) coordinate of bounding rectangle, in decimal degrees.
#' @param xmax Numeric. East longitude (x) coordinate of bounding rectangle, in decimal degrees.
#' @param ymin Numeric. South latitude (y) coordinate of bounding rectangle, in decimal degrees.
#' @param ymax Numeric. North latitude (y) coordinate of bounding rectangle, in decimal degrees.
#'
#' @return Returns an S4 object of class [rICESat2Veg::icesat2.atl08_dt]
#' containing the ATL08 Terrain attributes.
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
#atl08_h5<-ATL08_read(atl08_path=atl08_path)
#'
#'# Extracting ATL08-derived terrain Metrics
#'atl08_terrain_dt<-ATL08_terrain_attributes(atl08_h5=atl08_h5)
#'
#' # Bounding rectangle coordinates
#' xmin <- -107.7
#' xmax <- -106.5
#' ymin <- 32.75
#' ymax <- 42.75
#'
#'# Clipping ATL08-derived terrain metrics by boundary box extent
#'atl08_terrain_dt_clip <- ATL08_terrain_dt_clipBox(atl08_terrain_dt, xmin, xmax, ymin, ymax)
#'
#'close(atl08_h5)
#'@import hdf5r stats
#'@export
ATL08_terrain_dt_clipBox <- function(atl08_terrain_dt, xmin, xmax, ymin, ymax) {

  if (!class(atl08_terrain_dt)[1]=="icesat2.atl08_dt"){
    stop("atl08_terrain_dt needs to be an object of class 'icesat2.atl08_dt' ")
  }

  if (any(is.na(atl08_terrain_dt@dt))) {
    atl08_terrain_dt<-na.omit(atl08_terrain_dt@dt)
  } else {

    atl08_terrain_dt<-atl08_terrain_dt@dt
  }


  # xmin ymin xmax ymax
  mask <-
    atl08_terrain_dt$longitude >= xmin &
    atl08_terrain_dt$longitude <= xmax &
    atl08_terrain_dt$latitude >= ymin &
    atl08_terrain_dt$latitude <= ymax

  mask[!stats::complete.cases(mask)] <- FALSE
  mask <- (seq_along(atl08_terrain_dt$latitude))[mask]
  newFile <- atl08_terrain_dt[mask, ]

  newFile<- new("icesat2.atl08_dt", dt = newFile)

  # newFile<- new("gedi.level1b.dt", dt = level1bdt[mask,])
  return(newFile)
}
