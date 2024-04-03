#' Statistics of ATL03 and ATL08 joined photons attributes within a given area
#'
#' @description Computes a series of statistics ATL03 and ATL08 joined photons attributes within
#' area defined by a polygon
#'
#' @usage ATL03_ATL08_photons_attributes_dt_polyStat(atl03_atl08_dt, func, seg_length,ph_class,beam,quality_ph,night_flag)
#'
#' @param atl03_atl08_dt  An S4 object of class [ICESat2VegR::icesat2.atl08_dt] containing ATL03 and ATL08 data
#' (output of [ATL03_ATL08_photons_attributes_dt_join()] function).
#' @param func The function to be applied for computing the defined statistics
#' @param poly_id Polygon id. If defined, statistics will be computed for each polygon
#'
#' @return Returns an S4 object of class [ICESat2VegR::icesat2.atl08_dt]
#' Containing Statistics of ATL08 classified canopy photons
#'
#' @examples
#'
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
#'# Extracting ATL03 and ATL08 photons and heights
#'atl03_atl08_dt<-ATL03_ATL08_photons_attributes_dt_join(atl03_h5,atl08_h5)
#'head(atl03_atl08_dt)
#'
#'# Specifying the path to shapefile
#' polygon_filepath <- system.file("extdata", "polygon.shp", package = "ICESat2VegR")
#'
#'# Reading shapefile as sf object
#'polygon <- terra::vect(polygon_filepath)
#'
#'# Clipping ATL08 terrain attributes by Geometry
#' atl03_atl08_dt_clip <- ATL03_ATL08_joined_dt_clipGeometry(atl03_atl08_dt, polygon, split_by = "FID")
#'
#'# Computing the maximum ph_h by polygon id
#'max_ph_h <-ATL03_ATL08_photons_attributes_dt_polyStat(atl03_atl08_dt_clip, func=max(ph_h),poly_id="poly_id")
#'head(max_ph_h)
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
#' # Computing a series of ph_h statistics from customized function
#'ph_h_metrics <-ATL03_ATL08_photons_attributes_dt_polyStat(atl03_atl08_dt, func=mySetOfMetrics(ph_h),poly_id="poly_id")
#'head(ph_h_metrics)
#'
#'close(atl03_h5)
#'close(atl08_h5)
#' @import data.table lazyeval
#' @export
ATL03_ATL08_photons_attributes_dt_polyStat <- function(atl03_atl08_dt, func, poly_id=NULL) {

  if (!inherits(atl03_atl08_dt, "icesat2.atl03atl08_dt")){
    stop("atl03_atl08_dt needs to be an object of class 'icesat2.atl03atl08_dt' ")
  }

  if (any(is.na(atl03_atl08_dt))) {
    atl03_atl08_dt2<-na.omit(atl03_atl08_dt)
  } else {

    atl03_atl08_dt2<-atl03_atl08_dt
  }


  # Add data.table operator
  `:=` <- data.table::`:=`

  call <- substitute(func)

  if (is.null(poly_id)) {
    metrics <- lazy_apply_dt_call(dt = atl03_atl08_dt, call = call)
    metrics <- as.data.table(metrics)
    if (ncol(metrics) < 2) {
      colnames(metrics) <- paste0(call)[1]
    }
  } else {
    metrics <- lazy_apply_dt_call(dt = atl03_atl08_dt, call = call, group.by = paste0("by = ", poly_id))

    if (ncol(metrics) < 3) {
      colnames(metrics)[2] <- paste0(call)[1]
    }
  }

  return(metrics)
}
