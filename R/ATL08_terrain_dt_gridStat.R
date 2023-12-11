#' Statistics of ATL08 terrain attributes at grid level
#'
#' @description This function computes a series of user defined descriptive statistics within
#' each given grid cell for ATL08 terrain attributes
#'
#' @usage ATL08_terrain_segStat(atl03_atl08_dt, func, res)
#'
#' @param ATL08_terrain_dt  An S4 object of class [rICESat2Veg::icesat2.atl08_dt] containing ATL08 data
#' (output of [ATL08_terrain_attributes()] or [ATL08_terrain_segStat()] functions).
#' @param func The function to be applied for computing the defined statistics
#' @param res Spatial resolution in decimal degrees for the output SpatRast raster layer. Default is 0.5.
#'
#' @return Return a SpatRast raster layer(s) of selected ATL08 terrain attribute(s)
#'
#' @examples
#'
#'# Specifying the path to ATL08 file (zip file)
#'outdir = tempdir()
#'
#'atl08_zip <- system.file("extdata",
#'                   "ATL08_20220401221822_01501506_005_01.zip",
#'                   package="rICESat2Veg")
#'
#'# Unzipping ATL08 file
#'atl08_path <- unzip(atl08_zip,exdir = outdir)
#'
#'# Reading ATL08 data (h5 file)
#atl08_h5<-ATL08read(atl08_path=atl08_path)
#'
#'# Extracting ATL03 and ATL08 photons and heights
#'atl08_terrain_dt<-ATL08_terrain_attributes(atl08_h5)
#'head(atl08_terrain_dt)
#'
#'# Computing the mean h_te_best_fit at 30 m grid cell
#'mean_h_terrain <-ATL08_terrain_dt_gridStat(atl08_terrain_dt, func=mean(h_terrain), res=0.5)
#'
#'plot(top_terrain)
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
#' # Computing a series of h_te_best_fit statistics at 30 m grid cellfrom customized function
#'h_terrain_metrics <-ATL08_terrain_dt_gridStat(atl08_terrain_dt, func=mySetOfMetrics(h_te_best_fit),res=0.5)
#'
#'plot(h_terrain_metrics)
#'
#'close(atl03_h5)
#'close(atl08_h5)
#' @import data.table lazyeval
#' @export
ATL08_terrain_dt_gridStat <- function(atl08_terrain_dt, func, res = 0.5) {

  if (!class(atl08_terrain_dt)[1]=="icesat2.atl08_dt"){
    stop("ATL08_terrain_dt needs to be an object of class 'icesat2.atl08_dt' ")
  }

  if (any(is.na(atl08_terrain_dt@dt))) {
    atl08_terrain_dt2<-na.omit(atl08_terrain_dt@dt)
  } else {

    atl08_terrain_dt2<-atl08_terrain_dt@dt
  }

  if (!nrow(atl08_terrain_dt2)>1){
    stop(paste("ATL08_terrain_dt is invalid. It contain only", nrow(atl08_terrain_dt2),"observations"))
  }


  requireNamespace("data.table")
  cells <- NA

  # Add data.table operator
  `:=` <- data.table::`:=`

  call <- lazy_call(func)

  vect <- terra::vect(
    atl08_terrain_dt2,
    geom = c("longitude", "latitude"),
    crs = "epsg:4326"
  )
  layout <- terra::rast(terra::ext(vect), resolution = res, vals = NA, crs = "epsg:4326")

  atl08_terrain_dt2[, cells := terra::cells(layout, vect)[, 2]]
  metrics <- lazy_apply_dt_call(atl08_terrain_dt2, call, group.by = "by = cells")

  if (any(is.na(metrics))) {
  metrics<-na.omit(metrics)
  }

  n_metrics <- ncol(metrics) - 1
  bbox <- terra::ext(vect)
  output <-
    terra::rast(
      bbox,
      resolution = res,
      nlyrs = n_metrics,
      crs = "epsg:4326"
    )


  for (metric in 1:n_metrics) {
    output[[metric]][metrics$cells] <- metrics[[metric]]
  }

  return(output)
}

