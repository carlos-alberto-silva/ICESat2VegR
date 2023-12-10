#' Compute descriptive statistics of GEDI Elevation and Height Metrics
#'
#' @description Computes a Series of Statistics from GEDI derived Elevation and Height Metrics (Level2A)
#' within a given area defined or not by a polygon
#'
#' @usage polyStatsLevel2AM(level2AM, func, id=NULL)
#'
#' @param level2AM A GEDI Level2AM object (output of [getLevel2AM()] function).
#' An S4 object of class "data.table".
#' @param func The function to be applied for computing the defined statistics
#' @param id A vector containing the polygon id for each GEDI observation. Default is NULL
#'
#' @return Returns an S4 object of class [data.table::data.table]
#' Containing Statistics of GEDI level2A defined metrics
#'
#' @seealso \url{https://lpdaac.usgs.gov/products/gedi02_av002/}
#'
#' @examples
#' # Specifying the path to GEDI level2A data (zip file)
#' outdir <- tempdir()
#' level2A_fp_zip <- system.file("extdata",
#'   "GEDI02_A_2019108080338_O01964_T05337_02_001_01_sub.zip",
#'   package = "rGEDI"
#' )
#'
#' # Unzipping GEDI level2A data
#' level2Apath <- unzip(level2A_fp_zip, exdir = outdir)
#'
#' # Reading GEDI level2A data (h5 file)
#' level2a <- readLevel2A(level2Apath = level2Apath)
#'
#' # Specifying the path to shapefile
#' polygon_filepath <- system.file("extdata", "stands_cerrado.shp", package = "rGEDI")
#'
#' # Reading shapefile as sf object
#' library(sf)
#' polygon <- sf::st_read(polygon_filepath)
#'
#' # Extracting GEDI Eleveation and Relative Metrics (level2A)
#' level2AM <- getLevel2AM(level2a)
#' head(level2AM)
#'
#' # Clipping GEDI data by geometry
#' level2AM_clip <- clipLevel2AMGeometry(level2AM, polygon, split_by = "id")
#'
#' #' Define your own function
#' mySetOfMetrics <- function(x) {
#'   metrics <- list(
#'     min = min(x), # Min of x
#'     max = max(x), # Max of x
#'     mean = mean(x), # Mean of x
#'     sd = sd(x) # Sd of x
#'   )
#'   return(metrics)
#' }
#'
#' # Computing the maximum of RH100
#' RH100max <- polyStatsLevel2AM(level2AM_clip, func = max(rh100), id = NULL)
#'
#' # Computing the maximum of RH100 stratified by polygon
#' RH100max_poly <- polyStatsLevel2AM(level2AM_clip, func = max(rh100), id = NULL)
#'
#' # Computing a serie statistics for GEDI metrics stratified by polygon
#' RH100metrics <- polyStatsLevel2AM(level2AM_clip,
#'   func = mySetOfMetrics(rh100),
#'   id = level2AM_clip$id
#' )
#'
#' head(RH100metrics)
#'
#' close(level2a)
#' @import data.table lazyeval
#' @export
#'
#'

atl08_path<-"Z:\\01_Projects\\04_NASA_ICESat2\\10_others\\rICESat2Veg\\inst\\exdata\\ATL08_20220401221822_01501506_005_01.h5"
atl08_h5<-ATL08read(ATL08_path=atl08_path)
atl03_path<-"C:\\Users\\c.silva\\Documents\\rICESat2Veg\\inst\\exdata\\ATL03_20220401221822_01501506_005_01.h5"
atl03_h5<-ATL03read(atl03_path="C:\\Users\\c.silva\\Documents\\rICESat2Veg\\inst\\exdata\\ATL03_20220826063041_09981605_005_01.h5")


atl03_atl08_df<-ATL03_ATL08join(atl03_h5,atl08_h5, beam = "gt1l")

head(atl03_atl08_df)


mySetOfMetrics <- function(x) {
   metrics <- list(
     min = min(x), # Min of x
     max = max(x), # Max of x
     mean = mean(x), # Mean of x
     sd = sd(x) # Sd of x
   )
   return(metrics)
}
#'
#' # Computing the maximum of RH100
#' RH100max <- polyStatsLevel2AM(level2AM_clip, func = max(rh100), id = NULL)
#'

RH100max <-ATL08_canopy_segStat(atl03_atl08_df, func=max(ph_h), seg_lenghth = 30)


ATL08_canopy_segStat <- function(atl03_atl08_df, func, seg_lenghth = 30) {

  atl03_atl08_df2<-atl03_atl08_df[atl03_atl08_df$classed_pc_flag %in% c(1,2),]

  bin<-cut(atl03_atl08_df2$dist_ph_along, breaks=seq(0,max(atl03_atl08_df2$dist_ph_along),seg_lenghth),labels = FALSE)

  atl03_atl08_df2$bin<-bin

  requireNamespace("data.table")

  # Add data.table operator
  `:=` <- data.table::`:=`

  call <- lazy_call(func)

  if (is.null(bin)) {
    metrics <- lazy_apply_dt_call(dt = atl03_atl08_df2, call = call)
    metrics <- as.data.table(metrics)
    if (ncol(metrics) < 2) {
      colnames(metrics) <- paste0(call)[1]
    }
  } else {
    metrics <- lazy_apply_dt_call(dt = atl03_atl08_df2, call = call, group.by = "by = bin")
    if (ncol(metrics) < 3) {
      colnames(metrics)[2] <- paste0(call)[1]
    }
  }

  return(metrics)
}
