#' Statistics of ATL03 and ATL08 joined photons attributes within a given area
#'
#' @description Computes a series of statistics ATL03 and ATL08 joined photons attributes within
#' area defined by a polygon
#'
#' @usage ATL08_terrain_dt_polyStat(atl08_terrain_dt, func, poly_id)
#'
#' @param atl08_terrain_dt  An S4 object of class [rICESat2Veg::icesat2.atl08_dt] containing ATL08 terrain attributes
#' (output of [rICESat2Veg::ATL08_terrain_attribute_dt()] function).
#' @param func The function to be applied for computing the defined statistics
#' @param poly_id Polygon id. If defined, statistics will be computed for each polygon
#'
#' @return Returns an S4 object of class [rICESat2Veg::icesat2.atl08_dt]
#' Containing Statistics of ATL08 classified terrain photons
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
#atl08_h5<-ATL08read(atl08_path=atl08_path)
#'
#'# Extracting ATL08-derived terrain Metrics
#'atl08_terrain_dt<-ATL08_terrain_attributes(atl08_h5=atl08_h5)
#'
#' # Specifying the path to shapefile
#' polygon_filepath <- system.file("extdata", "polygon.shp", package = "rICESat2Veg")
#'
#' # Reading shapefile as sf object
#'polygon <- terra::vect(polygon_filepath)
#'
#' # Clipping ATL08 terrain attributes by Geometry
#' atl08_terrain_dt_clip <- ATL08_terrain_dt_clipGeometry(atl08_terrain_dt, polygon, split_by = "FID")
#'
#'# Computing the mean h_te_best_fit by polygon id
#'mean_h_te_best_fit <-ATL08_terrain_dt_polyStat(atl08_terrain_dt_clip, func=mean(h_te_best_fit),poly_id="poly_id")
#'head(mean_h_te_best_fit)
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
#'h_te_best_fit_metrics <-ATL08_terrain_dt_polyStat(atl08_dt, func=mySetOfMetrics(h_te_best_fit),poly_id="poly_id")
#'head(h_te_best_fit_metrics)
#'
#'close(atl03_h5)
#'close(atl08_h5)
#' @import data.table lazyeval
#' @export
ATL08_terrain_dt_polyStat <- function(atl08_terrain_dt, func, poly_id=NULL) {

  if (!class(atl08_terrain_dt)[1]=="icesat2.atl08_dt"){
    stop("atl08_terrain_dt needs to be an object of class 'icesat2.atl08_dt' ")
  }

  if (any(is.na(atl08_terrain_dt@dt))) {
    atl08_terrain_dt2<-na.omit(atl08_terrain_dt@dt)
  } else {

    atl08_terrain_dt2<-atl08_terrain_dt@dt
  }

  requireNamespace("data.table")

  # Add data.table operator
  `:=` <- data.table::`:=`

  call <- lazy_call(func)

  if (is.null(poly_id)) {
    metrics <- lazy_apply_dt_call(dt = atl08_terrain_dt@dt, call = call)
    metrics <- as.data.table(metrics)
    if (ncol(metrics) < 2) {
      colnames(metrics) <- paste0(call)[1]
    }
  } else {
    metrics <- lazy_apply_dt_call(dt = atl08_terrain_dt@dt, call = call, group.by = paste0("by = ", poly_id))

    if (ncol(metrics) < 3) {
      colnames(metrics)[2] <- paste0(call)[1]
    }
  }

  #metrics<- new("icesat2.atl08_dt", dt = metrics)
  return(metrics)
}
