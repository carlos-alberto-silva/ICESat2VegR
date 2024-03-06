randomSamplingWorker <- function(dt, size) {
  saneSample(dt, size)
}

spacedSamplingWorker <- function(dt, size, radius, spatialIndex = NULL) {
  if (is.null(spatialIndex)) {
    spatialIndex <- ANNIndex$new(dt$longitude, dt$latitude)
  }

  idx <- pkg_module$sampleMinDistanceRcpp(
    dt$longitude,
    dt$latitude,
    radius,
    size,
    spatialIndex$tree,
    as.integer(runif(1, 0, 2147483647))
  )

  dt[idx]
}
# grid_size = 0.1
# size = 0.999999
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
  .N <- data.table::.N
  # variable <- "h_canopy"

  y <- dt[, .SD, .SDcols = variable][[1]]
  hist_bins <- hist(y, plot = FALSE, ...)

  grouped <- dt[, cbind(.SD, breaks = cut(y, hist_bins$breaks))]
  grouped[, I := .I]

  grouped[, do.call(chainSampling$fn, c(list(.SD), chainSampling$params)), by = breaks]
}



#' @export
sample <- function(dt, size = NULL, method = NULL, ...) {
  allMethods <- gsub("sample\\.", "", as.character(.S3methods("sample")))
  if (any(class(dt) %in% allMethods)) {
    UseMethod("sample")
  }

  base::sample(dt, size, ...)
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

#' Sample spatial points given a minimum distance between them
#'
#' @param dt [`icesat2.atl08_dt-class`] data.table. Which can be extracted using
#' [`ATL08_seg_attributes_dt()`]
#' @param method the sampling method, either one of the sampling methods
#' provided [`randomSampling()`], [`spacedSampling()`], [`gridSampling()`],
#' [`stratifiedSampling()`]
#'
#' @return a subset sample of the input [`icesat2.atl08_dt-class`]
#'
#' @details
#' This function uses a wrap around C++ ANN library from Mount and Arya (2010).
#' for calculating the nearest neighbors and including them in a taboo list
#' which won't be considered for the next sampling.
#'
#' @seealso
#' Mount, D. M.; Arya, S. ANN: A Library for Approximate Nearest Neighbor Searching,
#' available in <https://www.cs.umd.edu/~mount/ANN/>
#'
#' @export
`sample.icesat2.atl08_dt` <- function(dt, method, ...) {
  do.call(method$fn, c(list(dt), method$params))
}




#' Sample spatial points given a minimum distance between them
#'
#' @param dt [`icesat2.atl03_seg_dt-class`] data.table. Which can be extracted using
#' [`ATL08_seg_attributes_dt()`]
#' @param sampling_method
#' between points. default NULL will calculate a new ANN index tree.
#'
#' @return a subset sample of the input [`icesat2.atl03_seg_dt-class`]
#'
#' @details
#' This function uses a wrap around C++ ANN library from Mount and Arya (2010).
#' for calculating the nearest neighbors and including them in a taboo list
#' which won't be considered for the next sampling.
#'
#' @seealso
#' Mount, D. M.; Arya, S. ANN: A Library for Approximate Nearest Neighbor Searching,
#' available in <https://www.cs.umd.edu/~mount/ANN/>
#'

#' @export
`sample.icesat2.atl03_seg_dt` <- function(dt, sampling_method, ...) {

}




########################
## SAMPLING METHODS
########################
setRefClass("icesat2_sampling_method")

#' @export
"+.icesat2_sampling_method" <- function(e1, e2) {
  e1$params[["chainSampling"]] <- e2
  e1
}

genericSamplingMethod <- function(fn, ..., chainSampling = NULL) {
  params <- list(...)
  params_string <- ""
  if (length(params) > 0) {
    params_names <- names(params)

    # Get params names whose values are not NULL
    params_with_value_names <- unlist(sapply(names(params), function(x) if (!is.null(params[[x]])) x))

    params_string <- paste(gsub("^", ", ", params_names), collapse = "")
    if (!is.null(params_with_value_names)) {
      params_string <- gsub(
        gettextf("(%s)", params_with_value_names),
        gettextf("\\1 = %s", list(params_with_value)),
        params_string
      )
    }
  }

  if (!is.null(chainSampling)) {
    params_string <- gettextf("%s, chainSampling = %s", params_string, list(chainSampling))
  }

  fn_definition <- sprintf('
  function(size%1$s, ...) {
    methodContainer <- list(
      fn = fn,
      params = list(size%1$s, ...)
    )

    prepend_class(methodContainer, "icesat2_sampling_method")
    methodContainer
  }
  ', params_string)

  params_with_value <-
    eval(parse(text = fn_definition))

  params_with_value
}

#' Pure random sampling method
#'
#' @param dt [`icesat2.atl03_seg_dt-class`] data.table. Which can be extracted using
#' [`ATL08_seg_attributes_dt()`]
#' @param size the sample size. Either an integer of absolute number of samples or
#' if it is between (0, 1) it will sample a percentage relative to the number of
#' available observations within the group.
#' @param ... generic params not used by randomSampling
#'
#' @return A [`icesat2_sampling_method-class`] for defining the method used in [`sample()`]
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
#'
#' @return A [`icesat2_sampling_method-class`] for defining the method used in [`sample()`]
#'
#' @export
gridSampling <- genericSamplingMethod(gridSamplingWorker, grid_size = NULL)



#' Get samples stratified by a variable binning histogram
#'
#' @param size the sample size. Either an integer of absolute number of samples or
#' if it is between (0, 1) it will sample a percentage relative to the number of
#' available observations within the group.
#' @param variable Variable name used for the stratification
#' @param chainSampling chains different methods of sampling by providing the result
#' of another samplingMethod [`randomSampling()`], [`spacedSampling()`], [`stratifiedSampling()`].
#' You can chain for example, a gridSampling(100, )
#' @param ...
#' forward to the [`graphics::hist()`], where you can manually define the breaks.
#'
#' @return A [`icesat2_sampling_method-class`] for defining the method used in [`sample()`]
#'
#' @export
stratifiedSampling <- genericSamplingMethod(stratifiedSamplingWorker, variable = NULL)

#' Get observations given a minimum radius distance between samples
#'
#' @param size the sample size. Either an integer of absolute number of samples or
#' if it is between (0, 1) it will sample a percentage relative to the number of
#' available observations within the group.
#' @param variable Variable name used for the stratification
#' @param ...
#' forward to the [`graphics::hist()`], where you can manually define the breaks.
#'
#' @return A [`icesat2_sampling_method-class`] for defining the method used in [`sample()`]
#'
#' @export
spacedSampling <- genericSamplingMethod(spacedSamplingWorker, radius = NULL)




# ggplot() +
#   geom_spatvector(data = geom, color = "orange", fill = NA) +
#   geom_point(data = sampled, aes(longitude, latitude, color = h_canopy), cex = 3) +
#   scale_x_continuous(breaks = seq(min(all_data$longitude), max(all_data$longitude), 1)) +
#   scale_y_continuous(breaks = seq(min(all_data$latitude), max(all_data$latitude), 1)) +
#   scale_color_gradient(low = "gray", high = "forestgreen") +
#   ggplot2::theme_bw()
