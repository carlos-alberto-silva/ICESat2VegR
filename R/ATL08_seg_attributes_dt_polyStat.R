#' Statistics of ATL03 and ATL08 joined photons attributes within a given area
#'
#' @description Computes a series of statistics ATL03 and ATL08 joined photons attributes within
#' area defined by a polygon
#'
#' @usage ATL08_seg_attributes_dt_polyStat(atl08_terrain_dt, func, poly_id)
#'
#' @param atl08_canopy_dt  An S4 object of class "icesat2.ATL08_canopy_dt" containing ATL08 canopy attributes
#' (output of [ATL08_canopy_attribute_dt()] function).
#' @param func The function to be applied for computing the defined statistics
#' @param poly_id Polygon id. If defined, statistics will be computed for each polygon
#'
#' @return Returns an S4 object of class [rICESat2Veg::icesat2.ATL08_canopy_dt]
#' Containing Statistics of ATL08 classified canopy photons
#'
#' @examples
#' # Specifying the path to ATL08 file (zip file)
#'outdir = tempdir()
#'atl08_zip <- system.file("extdata",
#'                   "ATL08_20220401221822_01501506_005_01.zip",
#'                   package="rICESat2Veg")
#'
#' # Unzipping ATL08 file
#'atl08_path <- unzip(atl08_zip,exdir = outdir)
#'
#' # Reading ATL08 data (h5 file)
#atl08_h5<-ATL08_read(ATL08_path=atl08_path)
#'
#'# Extracting ATL08-derived canopy Metrics
#'atl08_canopy_dt<-ATL08_seg_attributes_dt(atl08_h5=atl08_h5)
#'
#' # Specifying the path to shapefile
#' polygon_filepath <- system.file("extdata", "polygon.shp", package = "rICESat2Veg")
#'
#' # Reading shapefile as sf object
#'polygon <- terra::vect(polygon_filepath)
#'
#' # Clipping ATL08 canopy attributes by Geometry
#' atl08_canopy_dt_clip <- ATL08_seg_attributes_dt_clipGeometry(atl08_canopy_dt, polygon, split_by = "FID")
#'
#'# Computing the max h_canopy by polygon id
#'max_h_canopy <-ATL08_seg_attributes_dt_polyStat(atl08_canopy_dt_clip, func=max(h_canopy),poly_id="poly_id")
#'head(max_h_canopy)
#'
#'# Define your own function
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
#' # Computing a series of canopy statistics from customized function
#'h_canopy_metrics <-ATL08_seg_attributes_dt_polyStat(ATL08_canopy_dt, func=mySetOfMetrics(h_canopy),poly_id="poly_id")
#'head(h_canopy_metrics)
#'
#'close(atl03_h5)
#'close(atl08_h5)
#' @import data.table lazyeval
#' @export
ATL08_seg_attributes_dt_polyStat <- function(atl08_canopy_dt, func, poly_id=NULL) {

  if (!class(atl08_canopy_dt)[1]=="icesat2.ATL08_canopy_dt"){
    stop("atl08_canopy_dt needs to be an object of class 'icesat2.ATL08_canopy_dt' ")
  }

  if (any(is.na(atl08_canopy_dt@dt))) {
    atl08_canopy_dt2<-na.omit(atl08_canopy_dt@dt)
  } else {

    atl08_canopy_dt2<-atl08_canopy_dt@dt
  }

  
  # Add data.table operator
  `:=` <- data.table::`:=`

  call <- lazy_call(func)

  if (is.null(poly_id)) {
    metrics <- lazy_apply_dt_call(dt = atl08_canopy_dt@dt, call = call)
    metrics <- as.data.table(metrics)
    if (ncol(metrics) < 2) {
      colnames(metrics) <- paste0(call)[1]
    }
  } else {
    metrics <- lazy_apply_dt_call(dt = atl08_canopy_dt@dt, call = call, group.by = paste0("by = ", poly_id))

    if (ncol(metrics) < 3) {
      colnames(metrics)[2] <- paste0(call)[1]
    }
  }

  #metrics<- new("icesat2.ATL08_canopy_dt", dt = metrics)
  return(metrics)
}
