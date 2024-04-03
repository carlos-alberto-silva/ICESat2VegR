#' Statistics of ATL08 Terrain and Canopy Attributes at grid level
#'
#' @description This function computes a series of user defined descriptive statistics within
#' each grid cell for ATL08 Terrain and Canopy Attributes
#'
#' @usage ATL08_seg_attributes_dt_gridStat(atl08_seg_att_dt, func, res)
#'
#' @param atl08_seg_att_dt  An S4 object of class [ICESat2VegR::icesat2.atl08_dt] containing ATL08 data
#' (output of [ATL08_seg_attributes_dt()] functions).
#' @param func The function to be applied for computing the defined statistics
#' @param res Spatial resolution in decimal degrees for the output SpatRast raster layer. Default is 0.5.
#'
#' @return Return a SpatRast raster layer(s) of selected ATL08 terrain and canopy attribute(s)
#'
#' @examples
#' # Specifying the path to ATL08 file
#' atl08_path <- system.file("extdata",
#'   "atl08_clip.h5",
#'   package = "ICESat2VegR"
#' )
#'
#' # Reading ATL08 data (h5 file)
#' atl08_h5 <- ATL08_read(atl08_path = atl08_path)
#'
#' # Extracting ATL08-derived terrain and canopy attributes
#' atl08_seg_att_dt <- ATL08_seg_attributes_dt(atl08_h5 = atl08_h5)
#'
#' # Computing the top h_canopy at 0.05 degree grid cell
#' res <- 0.0001
#' mean_h_canopy <- ATL08_seg_attributes_dt_gridStat(atl08_seg_att_dt, func = mean(h_canopy), res = res)
#'
#' plot(mean_h_canopy)
#'
#' # Define your own function
#' mySetOfMetrics <- function(x) {
#'   metrics <- list(
#'     min = min(x, na.rm = TRUE), # Min of x
#'     max = max(x, na.rm = TRUE), # Max of x
#'     mean = mean(x, na.rm = TRUE), # Mean of x
#'     sd = sd(x, na.rm = TRUE) # Sd of x
#'   )
#'   return(metrics)
#' }
#'
#' res <- 0.05
#' # Computing a series of h_canopy statistics at 0.05 degree grid cell from customized function
#' h_canopy_metrics <- ATL08_seg_attributes_dt_gridStat(atl08_seg_att_dt, func = ~ mySetOfMetrics(h_canopy), res = res)
#'
#' plot(h_canopy_metrics)
#'
#' close(atl08_h5)
#' @import data.table lazyeval
#' @export
ATL08_seg_attributes_dt_gridStat <- function(atl08_seg_att_dt, func, res = 0.5) {
  if (!inherits(atl08_seg_att_dt, "icesat2.atl08_dt")) {
    stop("atl08_seg_att_dt needs to be an object of class 'icesat2.atl08_dt' ")
  }

  if (any(is.na(atl08_seg_att_dt))) {
    atl08_seg_att_dt2 <- na.omit(atl08_seg_att_dt)
  } else {
    atl08_seg_att_dt2 <- atl08_seg_att_dt
  }

  if (!nrow(atl08_seg_att_dt2) > 1) {
    stop(paste("atl08_seg_att_dt is invalid. It contain only", nrow(atl08_seg_att_dt2), "observations"))
  }


  cells <- NA

  # Add data.table operator
  `:=` <- data.table::`:=`

  call <- substitute(func)

  vect <- to_vect(atl08_seg_att_dt2)
  layout <- terra::rast(terra::ext(vect), resolution = res, vals = NA, crs = "epsg:4326")

  suppressWarnings(atl08_seg_att_dt2[, cells := terra::cells(layout, vect)[, 2]])
  metrics <- lazy_apply_dt_call(atl08_seg_att_dt2, call, group.by = "by = cells")

  n_metrics <- ncol(metrics) - 1
  bbox <- terra::ext(vect)
  output <-
    terra::rast(
      bbox,
      resolution = res,
      nlyrs = n_metrics,
      crs = "epsg:4326"
    )


  names(output) <- names(metrics)[-1]
  metrics <- metrics[!is.nan(cells)]
  
  for (metric in 1:n_metrics) {
    output[[metric]][metrics$cells] <- metrics[[metric + 1]]
  }

  return(output)
}
