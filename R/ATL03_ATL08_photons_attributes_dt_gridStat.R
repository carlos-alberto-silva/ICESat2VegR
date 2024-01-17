#' Statistics of ATL03 and ATL08 photon attributes
#'
#' @description This function computes a series of user defined descriptive statistics within
#' each given grid cell for ATL03 and ATL08 photon attributes
#'
#' @usage ATL03_ATL08_photons_attributes_dt_gridStat(atl03_atl08_dt, func, res)
#'
#' @param atl03_atl08_dt  An S4 object of class [ICESat2VegR::icesat2.atl03_atl08_dt] containing ATL03 and ATL08  attributes
#' (output of the [ATL03_ATL08_photons_attributes_dt_join()] function).
#' @param func The function to be applied for computing the defined statistics
#' @param res Spatial resolution in decimal degrees for the output SpatRast raster layer. Default is 0.5.
#' @param ph_class Character vector indicating photons to process based on the classification (1=ground, 2=canopy, 3=top canopy),
#' Default is c(2,3)
#' @param beam Character vector indicating beams to process. Default is c("gt1l", "gt1r", "gt2l", "gt2r", "gt3l", "gt3r")
#' @param quality_ph Indicates the quality of the associated photon. 0=nominal, 1=possible_afterpulse, 2=possible_impulse_response_
#' effect, 3=possible_tep. Default is 0
#' @param night_flag Flag indicating the data were acquired in night conditions: 0=day, 1=night. Default is 1
#'
#' @return Return a SpatRast raster layer(s) of selected ATL03 and ATL08 photon attribute(s)
#'
#' @examples
#' outdir <- tempdir()
#' atl03_zip <- system.file("extdata",
#'   "ATL03_20220401221822_01501506_005_01.zip",
#'   package = "ICESat2VegR"
#' )
#'
#' atl08_zip <- system.file("extdata",
#'   "ATL08_20220401221822_01501506_005_01.zip",
#'   package = "ICESat2VegR"
#' )
#'
#' # Unzipping ATL03 file
#' atl03_path <- unzip(atl03_zip, exdir = outdir)
#'
#' # Unzipping ATL08 file
#' atl08_path <- unzip(atl08_zip, exdir = outdir)

#' # Reading ATL03 data (h5 file)
#' atl03_h5 <- ATL08_read(atl03_path = atl03_path)
#'
#' # Reading ATL08 data (h5 file)
#' atl08_h5 <- ATL08_read(atl08_path = atl08_path)
#'
#' # # Extracting ATL03 and ATL08 photons and heights
#' atl03_atl08_dt <- ATL03_ATL08_photons_attributes_dt_join(atl03_h5, atl08_h5)
#' head(atl03_atl08_dt)
#'
#' # Computing the mean of ph_h attribute at 0.05 degree grid cell
#' mean_ph_h <- ATL03_ATL08_photons_attributes_dt_gridStat(atl03_atl08_dt, func = mean(ph_h), res = 0.005)
#'
#' plot(mean_ph_h, xlim=c(-107.5,-106.5),ylim=c(38,39))
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
#' # Computing a series of ph_h stats at 0.05 degree grid cell from customized function
#' ph_h_metrics <- ATL03_ATL08_photons_attributes_dt_gridStat(atl03_atl08_dt, func = mySetOfMetrics(ph_h), res = 0.005)
#'
#' plot(ph_h_metrics, xlim=c(-107.5,-106.5),ylim=c(38,39))
#'
#' close(atl03_h5)
#' close(atl08_h5)
#' @import data.table lazyeval
#' @export
ATL03_ATL08_photons_attributes_dt_gridStat <- function(atl03_atl08_dt,
                                           func,
                                           res = 0.5,
                                           ph_class = c(2, 3),
                                           beam = c("gt1l", "gt1r", "gt2l", "gt2r", "gt3l", "gt3r"),
                                           quality_ph = 0,
                                           night_flag = 1) {
  if (!class(atl03_atl08_dt)[1] == "icesat2.atl03atl08_dt") {
    stop("atl03_atl08_dt needs to be an object of class 'icesat2.atl03atl08_dt' ")
  }

  if (any(is.na(atl03_atl08_dt))) {
    atl03_atl08_dt2 <- na.omit(atl03_atl08_dt)
  } else {
    atl03_atl08_dt2 <- atl03_atl08_dt
  }

  atl03_atl08_dt2 <- atl03_atl08_dt2[
    atl03_atl08_dt2$classed_pc_flag %in% ph_class &
      atl03_atl08_dt2$quality_ph %in% quality_ph &
      atl03_atl08_dt2$beam %in% beam &
      atl03_atl08_dt2$night_flag %in% night_flag,
  ]

  if (!nrow(atl03_atl08_dt2) > 1) {
    stop(paste("atl03_atl08_dt is invalid. It contain only", nrow(atl03_atl08_dt2), "observations"))
  }

  cells <- NA

  # Add data.table operator
  `:=` <- data.table::`:=`

  call <- substitute(func)

  vect <- terra::vect(
    atl03_atl08_dt2,
    geom = c("lon_ph", "lat_ph"),
    crs = "epsg:4326"
  )
  layout <- terra::rast(terra::ext(vect), resolution = res, vals = NA, crs = "epsg:4326")

  atl03_atl08_dt2[, cells := terra::cells(layout, vect)[, 2]]
  metrics <- atl03_atl08_dt2[, eval(call), by = cells]

  if (any(is.na(metrics))) {
    metrics <- na.omit(metrics)
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
    output[[metric]][metrics$cells] <- metrics[[metric+1]]
  }

  names(output)<-paste0(paste0(call)[2],"_",names(metrics)[-1])
  return(output)
}
