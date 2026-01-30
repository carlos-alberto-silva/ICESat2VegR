#' Sample method for applying multiple sampling methods
#'
#' @param x the generic input data to be sampled
#' @param ... generic params to pass ahead to the specific sampling function
#' @param method the sampling method to use.
#'
#' @details
#' It is expected that the user pass a `method` parameter within `...`
#'
#' @seealso
#' [`randomSampling()`], [`spacedSampling()`], [`gridSampling()`],
#' [`stratifiedSampling()`], [`geomSampling()`], [`rasterSampling()`]
#' @export
sample <- function(x, ..., method) {
  UseMethod("sample")
}

#' @keywords internal
sample.default <- function(x, ...) base::sample(x, ...)

# #' Sampling function which accepts a method for sampling
# #'
# #' @param x [`ICESat2VegR::icesat2.atl08_dt-class`] data.table. Which can be extracted using
# #' [`ATL08_seg_attributes_dt()`]
# #' @param method the sampling method, either one of the sampling methods
# #' provided [`randomSampling()`], [`spacedSampling()`], [`gridSampling()`],
# #' [`stratifiedSampling()`]
# #'
# #' @return a subset sample of the input [`ICESat2VegR::icesat2.atl08_dt-class`] plus addiotional
# #' columns regarding the grouping
# #'
# #' @seealso
# #' [`randomSampling()`], [`spacedSampling()`], [`gridSampling()`],
# #' [`stratifiedSampling()`], [`geomSampling()`], [`rasterSampling()`]
# #'
#' @keywords internal
`sample.icesat2.atl08_dt` <- function(x, ..., method) {
  do.call(method$fn, c(list(x), method$params))
}

# #' Sampling function which accepts a method for sampling
# #'
# #' @param x [`ICESat2VegR::icesat2.atl03_seg_dt-class`] data.table. Which can be extracted using
# #' [`ATL08_seg_attributes_dt()`]
# #' @param ... Not used, only for compatibility with generic method.
# #' @param method The sampling method, either one or a combination of the sampling methods
# #' provided [`randomSampling()`], [`spacedSampling()`], [`gridSampling()`],
# #' [`stratifiedSampling()`], [`geomSampling()`], [`rasterSampling()`]
# #'
# #' @return a subset sample of the input [`ICESat2VegR::icesat2.atl03_seg_dt-class`]
# #'
# #' @seealso
# #' [`randomSampling()`], [`spacedSampling()`], [`gridSampling()`],
# #' [`stratifiedSampling()`], [`geomSampling()`], [`rasterSampling()`]
# #'
#' @keywords internal
`sample.icesat2.atl03_seg_dt` <- function(x, ..., method) {
  method <- list(...)[["method"]]
  do.call(method$fn, c(list(x), method$params))
}

#' @keywords internal
`sample.icesat2.atl03_atl08_seg_dt` <- function(x, ..., method) {
  do.call(method$fn, c(list(x), method$params))
}

#' Class for sampling methods to be passed on for the sample function
#' @keywords internal
setRefClass("icesat2_sampling_method")

#' @keywords internal
"+.icesat2_sampling_method" <- function(e1, e2) {
  ii <- 1
  chainSampling <- e1
  chainSamplingDeepString <- list("$params$chainSampling")

  while (!is.null(chainSampling$params$chainSampling)) {
    chainSampling <- chainSampling$params$chainSampling
    ii <- ii + 1
    chainSamplingDeepString[""] <- "$params$chainSampling"
  }
  finalChainingString <- paste(chainSamplingDeepString, collapse = "")
  evaluater <- paste0("e1", finalChainingString, " <- e2", collapse = "")
  eval(parse(text = evaluater))
  return(invisible(e1))
}


#################################################################
## SAMPLING METHODS
#################################################################
# WORKERS
#
# Workers are the functions that actually do the sampling work.
#
# The sampling functions will just wrap around the worker so
# that we can chain them before actually doing the computation.
#################################################################
randomSamplingWorker <- function(dt, size) {
  saneSample(dt, size)
}

spacedSamplingWorker <- function(dt, size, radius, spatialIndex = NULL, chainSampling = NULL) {
  if (is.null(spatialIndex)) {
    spatialIndex <- ANNIndex$new(dt$longitude, dt$latitude)
  }

  if (size < 1 && size > 0) {
    size <- as.integer(round(nrow(dt) * size))
  }

  idx <- pkg_module$sampleMinDistanceRcpp(
    dt$longitude,
    dt$latitude,
    radius,
    size,
    spatialIndex$tree,
    as.integer(runif(1, 0, 2147483647))
  )

  if (is.null(chainSampling)) {
    return(dt[idx,])
  }
  return(do.call(chainSampling$fn, c(list(dt[idx,]), chainSampling$params)))
}

gridSamplingWorker <- function(dt, size, grid_size, chainSampling = NULL) {
  if (is.null(chainSampling)) chainSampling <- randomSampling(size)
  longitude <- latitude <- group <- NA
  x_grid <- dt[, as.integer((longitude - min(longitude)) / grid_size)]
  y_grid <- dt[, as.integer((latitude - min(latitude)) / grid_size)]

  grouped <- dt[, cbind(.SD, x_grid, y_grid)]
  grouped[, I := .I]
  grouped[, do.call(chainSampling$fn, c(list(.SD), chainSampling$params)), by = list(x_grid, y_grid)]
}

stratifiedSamplingWorker <- function(dt, size, variable, chainSampling = NULL, ...) {
  if (is.null(chainSampling)) chainSampling <- randomSampling(size)
  breaks <- NA
  .SD <- data.table::.SD
  # variable <- "h_canopy"

  y <- dt[, .SD, .SDcols = variable][[1]]
  hist_bins <- graphics::hist(y, plot = FALSE, ...)

  grouped <- dt[, cbind(.SD, breaks = cut(y, hist_bins$breaks))]
  grouped[, I := .I]

  grouped[, do.call(chainSampling$fn, c(list(.SD), chainSampling$params)), by = breaks]
}

geomSamplingWorker <- function(dt, size, geom, split_id = NULL, chainSampling = NULL) {
  geom_group <- NA
  `:=` <- data.table::`:=`
  .SD <- data.table::.SD

  if (is.null(chainSampling)) chainSampling <- randomSampling(size)

  if (is.null(split_id)) {
    split_id <- "rowid"
    geom[["rowid"]] <- seq_along(geom)
    geom[[names(geom)]]
  }

  v <- terra::vect(
    as.data.frame(dt),
    geom = c("longitude", "latitude"),
    crs = "epsg:4326"
  )

  geom_v <- terra::intersect(v, geom)
  dt[, geom_group := geom_v[[split_id]]]
  dt[, do.call(chainSampling$fn, c(list(.SD), chainSampling$params)), by = geom_group]
}

rasterSamplingWorker <- function(dt, size, raster, chainSampling = NULL) {
  raster_group <- NA

  if (is.null(chainSampling)) chainSampling <- randomSampling(size)

  v <- terra::vect(
    as.data.frame(dt),
    geom = c("longitude", "latitude"),
    crs = "epsg:4326"
  )

  v_rast <- terra::extract(raster, v)

  dt[, raster_group := as.factor(v_rast[, 2])]
  dt <- dt[!is.na(raster_group)]

  dt[, do.call(chainSampling$fn, c(list(.SD), chainSampling$params)), by = raster_group]
}

saneSample <- function(dt, n) {
  nrows <- nrow(dt)
  if (n < 1 && n > 0) {
    return(dt[sample(dt[, .I], round(n * nrows))])
  }
  if (nrows < n) {
    warning(sprintf(
      "Not enough samples within group, truncated to %d samples",
      nrows
    ))
    n <- nrows
  }
  return(dt[sample(dt[, .I], n)])
}



#################################################################
# Generic dynamic evaluated function
#
# This just creates the appropriate function to wrap around the
# workers based on their declared parameters
##################################################################
genericSamplingMethod <- function(fn) {
  fn_name <- as.character(substitute(fn))
  # formals extract params of a function as list
  params <- formals(fn)
  params[[1]] <- NULL
  params_string <- paste0(sapply(names(params), function(x) {
    if (!is.null(params[[x]]) && params[[x]] == "") {
      x
    } else {
      gettextf(
        "%s = %s", x,
        list(params[[x]])
      )
    }
  }), collapse = ", ")

  fn_definition <- sprintf('
  function(%1$s) {
    methodContainer <- list(
      fn = fn,
      params = list(%1$s),
      methodName = "%2$s"
    )

    prepend_class(methodContainer, "icesat2_sampling_method")
    methodContainer
  }
  ', params_string, fn_name)

  params_with_value <-
    eval(parse(text = fn_definition))

  params_with_value
}

##########################################################
# SAMPLING FUNCTIONS using the generic function above
##########################################################
#' Pure random sampling method
#'
#' @param size the sample size. Either an integer of absolute number of samples or
#' if it is between (0, 1) it will sample a percentage relative to the number of
#' available observations within the group.
#'
#' @return A [`ICESat2VegR::icesat2_sampling_method-class`] for defining the method used in [`sample()`]
#'
#' @export
randomSampling <- genericSamplingMethod(randomSamplingWorker)

#' Get samples stratified by grid cells of specified size
#'
#' @param size the sample size. Either an integer of absolute number of samples or
#' if it is between (0, 1) it will sample a percentage relative to the number of
#' available observations within the group.
#' @param grid_size generic params. The [`gridSampling()`] will take a `grid_size` parameter
#' defining the grid size for the sampling
#' @param chainSampling chains different methods of sampling by providing the result
#' of another samplingMethod [`randomSampling()`], [`spacedSampling()`], [`stratifiedSampling()`].

#'
#' @return A [`ICESat2VegR::icesat2_sampling_method-class`] for defining the method used in [`sample()`]
#'
#' @export
gridSampling <- genericSamplingMethod(gridSamplingWorker)



#' Get samples stratified by a variable binning histogram
#'
#' @param size the sample size. Either an integer of absolute number of samples or
#' if it is between (0, 1) it will sample a percentage relative to the number of
#' available observations within the group.
#' @param variable Variable name used for the stratification
#' @param chainSampling chains different methods of sampling by providing the result
#' of another samplingMethod [`randomSampling()`], [`spacedSampling()`], [`stratifiedSampling()`].
#' @param ...
#' forward to the [`graphics::hist()`], where you can manually define the breaks.
#'
#' @return A [`ICESat2VegR::icesat2_sampling_method-class`] for defining the method used in [`sample()`]
#'
#' @export
stratifiedSampling <- genericSamplingMethod(stratifiedSamplingWorker)

#' Get observations given a minimum radius distance between samples
#'
#' @param size the sample size. Either an integer of absolute number of samples or
#' if it is between (0, 1) it will sample a percentage relative to the number of
#' available observations within the group.
#' @param radius the minimum radius between samples.
#' @param spatialIndex optional parameter. You can create a spatial index for accelerating
#' the search space with [`ANNIndex`] and reuse that if needed elsewhere. It will create
#' one automatically everytime if you don't specify one.
#' @param chainSampling chains different methods of sampling by providing the result
#' of another samplingMethod [`randomSampling()`], [`spacedSampling()`], [`stratifiedSampling()`].
#'
#' @return A [`ICESat2VegR::icesat2_sampling_method-class`] for defining the method used in [`sample()`]
#'
#' @include ANNIndex.R
#' @export
spacedSampling <- genericSamplingMethod(spacedSamplingWorker)

#' Get observations given a minimum radius distance between samples
#'
#' @param size the sample size. Either an integer of absolute number of samples or
#' if it is between (0, 1) it will sample a percentage relative to the number of
#' available observations within the group.
#' @param geom a [`terra::SpatVector-class`] object opened with [`terra::vect()`]
#' @param split_id character. The variable name of the geometry to use as split factor
#' for sampling
#' @param chainSampling chains different methods of sampling by providing the result
#' of another samplingMethod [`randomSampling()`], [`spacedSampling()`], [`stratifiedSampling()`].
#'
#' @return A [`ICESat2VegR::icesat2_sampling_method-class`] for defining the method used in [`sample()`]
#'
#' @export
geomSampling <- genericSamplingMethod(geomSamplingWorker)

#' Get observations given a minimum radius distance between samples
#'
#' @param size the sample size. Either an integer of absolute number of samples or
#' if it is between (0, 1) it will sample a percentage relative to the number of
#' available observations within the group.
#' @param raster a [`terra::SpatRaster-class`] object opened with [`terra::rast()`]
#' @param chainSampling chains different methods of sampling by providing the result
#' of another samplingMethod [`randomSampling()`], [`spacedSampling()`], [`stratifiedSampling()`].
#'
#' @return A [`ICESat2VegR::icesat2_sampling_method-class`] for defining the method used in [`sample()`]
#'
#' @export
rasterSampling <- genericSamplingMethod(rasterSamplingWorker)


# ggplot() +
#   geom_spatvector(data = geom, color = "orange", fill = NA) +
#   geom_point(data = sampled, aes(longitude, latitude, color = h_canopy), cex = 3) +
#   scale_x_continuous(breaks = seq(min(all_data$longitude), max(all_data$longitude), 1)) +
#   scale_y_continuous(breaks = seq(min(all_data$latitude), max(all_data$latitude), 1)) +
#   scale_color_gradient(low = "gray", high = "forestgreen") +
#   ggplot2::theme_bw()
