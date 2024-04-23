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
#'
#' # Unzipping ATL03 file
#' atl03_path <- system.file("extdata",
#'   "atl03_clip.h5",
#'   package = "ICESat2VegR"
#' )
#'
#' # Unzipping ATL08 file
#' atl08_path <- system.file("extdata",
#'   "atl08_clip.h5",
#'   package = "ICESat2VegR"
#' )
#' # Reading ATL03 data (h5 file)
#' atl03_h5 <- ATL03_read(atl03_path = atl03_path)
#'
#' # Reading ATL08 data (h5 file)
#' atl08_h5 <- ATL08_read(atl08_path = atl08_path)
#'
#' # # Extracting ATL03 and ATL08 photons and heights
#' atl03_atl08_dt <- ATL03_ATL08_photons_attributes_dt_join(atl03_h5, atl08_h5)
#'
#' # Computing the mean of ph_h attribute at 0.0002 degree grid cell
#' mean_ph_h <- ATL03_ATL08_photons_attributes_dt_gridStat(atl03_atl08_dt, func = mean(ph_h), res = 0.0002)
#'
#' plot(mean_ph_h)
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
#' # Computing a series of ph_h stats at 0.0002 degree grid cell from customized function
#' ph_h_metrics <- ATL03_ATL08_photons_attributes_dt_gridStat(atl03_atl08_dt, func = mySetOfMetrics(ph_h), res = 0.0002)
#'
#' plot(ph_h_metrics)
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
  if (!inherits(atl03_atl08_dt, "icesat2.atl03atl08_dt")) {
    stop("atl03_atl08_dt needs to be an object of class 'icesat2.atl03atl08_dt' ")
  }

  atl03_atl08_dt2 <- na.omit(atl03_atl08_dt)

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

  vect <- to_vect(atl03_atl08_dt2)
  layout <- terra::rast(terra::ext(vect) + res/4, resolution = res, vals = NA, crs = "epsg:4326")

  atl03_atl08_dt2[, cells := terra::cells(layout, vect)[, 2]]
  metrics <- atl03_atl08_dt2[, eval(call), by = cells]

  if (any(is.na(metrics))) {
    metrics <- na.omit(metrics)
  }

  n_metrics <- ncol(metrics) - 1
  bbox <- terra::ext(vect) + res/4
  output <-
    terra::rast(
      bbox,
      resolution = res,
      nlyrs = n_metrics,
      crs = "epsg:4326"
    )


  for (metric in 1:n_metrics) {
    output[[metric]][metrics$cells] <- metrics[[metric + 1]]
  }

  names(output) <- names(metrics)[-1]
  return(output)
}
