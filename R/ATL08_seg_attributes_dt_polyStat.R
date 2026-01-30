#' Statistics of ATL08 Terrain and Canopy Attributes by Geometry
#'
#' @description Computes a series of statistics of ATL08 terrain and canopy attributes within
#' area defined by a polygon
#'
#' @param atl08_seg_att_dt  An S4 object of class
#' [`ICESat2VegR::icesat2.atl08_dt-class`] containing
#' ATL08 terrain and canopy attributes (output of
#' [`ICESat2VegR::ATL08_seg_attributes_dt()`] function).
#' @param func The function to be applied for computing the defined statistics
#' @param poly_id Polygon id. If defined, statistics will be computed for each polygon
#'
#' @return Returns an S4 object of class [`ICESat2VegR::icesat2.atl08_dt-class`]
#' Containing Statistics of ATL08 terrain and canopy attributes
#'
#' @examples
#' # Specifying the path to ATL08 file
#' atl08_path <- system.file(
#'   "extdata",
#'   "atl08_clip.h5",
#'   package = "ICESat2VegR"
#' )
#'
#' # Reading ATL08 data (h5 file)
#' atl08_h5 <- ATL08_read(atl08_path = atl08_path)
#'
#' # Extracting ATL08 terrain and canopy attributes
#' atl08_seg_att_dt <- ATL08_seg_attributes_dt(atl08_h5 = atl08_h5)
#'
#' # Specifying the path to shapefile
#' polygon_filepath <- system.file(
#'   "extdata",
#'   "clip_geom.shp",
#'   package = "ICESat2VegR"
#' )
#'
#' # Reading shapefile as sf object
#' polygon <- terra::vect(polygon_filepath)
#'
#' # Clipping ATL08 terrain and canopy attributes by Geometry
#' atl08_seg_att_dt_clip <- ATL08_seg_attributes_dt_clipGeometry(
#'   atl08_seg_att_dt,
#'   polygon,
#'   split_by = "id"
#' )
#'
#' # Computing the max h_canopy by polygon id
#' max_h_canopy <- ATL08_seg_attributes_dt_polyStat(
#'   atl08_seg_att_dt_clip,
#'   func = max(h_canopy),
#'   poly_id = "poly_id"
#' )
#' head(max_h_canopy)
#'
#' # Define your own function
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
#' h_canopy_metrics <- ATL08_seg_attributes_dt_polyStat(
#'   atl08_seg_att_dt_clip,
#'   func = mySetOfMetrics(h_canopy),
#'   poly_id = "poly_id"
#' )
#'
#' head(h_canopy_metrics)
#'
#' close(atl08_h5)
#' @import data.table
#' @export
ATL08_seg_attributes_dt_polyStat <- function(atl08_seg_att_dt, func, poly_id = NULL) {
  if (!inherits(atl08_seg_att_dt, "icesat2.atl08_dt")) {
    stop("atl08_seg_att_dt needs to be an object of class 'icesat2.atl08_dt' ")
  }

  # Add data.table operator
  `:=` <- data.table::`:=`

  call <- substitute(func)

  if (is.null(poly_id)) {
    metrics <- lazy_apply_dt_call(dt = atl08_seg_att_dt, call = call)
    metrics <- as.data.table(metrics)
    if (ncol(metrics) < 2) {
      colnames(metrics) <- paste0(call)[1]
    }
  } else {
    metrics <- lazy_apply_dt_call(dt = atl08_seg_att_dt, call = call, group.by = paste0("by = ", poly_id))

    if (ncol(metrics) < 3) {
      colnames(metrics)[2] <- paste0(call)[1]
    }
  }

  prepend_class(metrics, c("icesat2.atl08_dt","data.table","data.frame"))
  return(metrics)
}
