randomSamplingWorker <- function(dt, size) {
  dt[saneSample(seq_len(len = nrow(all_data)), size)]
}

radiusSamplingWorker <- function(dt, size, radius, spatialIndex = NULL) {
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
# all_data$nid <- NULL
# dt <- all_data
gridSamplingWorker <- function(dt, size, grid_size) {
  longitude <- latitude <- group <- NA
  x_grid <- dt[, as.integer((longitude - min(longitude)) / grid_size)]
  y_grid <- dt[, as.integer((latitude - min(latitude)) / grid_size)]

  grouped <- dt[, cbind(.SD, x_grid, y_grid)]
  grouped[, I := .I]
  idx <- grouped[, list(idx = saneSample(I, size)), by = list(x_grid, y_grid)]

  dt[idx$idx, cbind(.SD, idx[, .SD, .SDcols = -"idx"])]
}


stratifiedSamplingWorker <- function(dt, size, variable, ...) {
  group <- NA
  .N <- data.table::.N
  # variable <- "h_canopy"

  hist_bins <- hist(dt[, get(variable)], plot = FALSE, ...)

  grouped <- dt[, cbind(.SD, breaks = cut(get(variable), hist_bins$breaks))]
  grouped[, I := .I]

  idx <- grouped[, list(idx = saneSample(I, size)), by = breaks]
  dt[idx$idx, cbind(.SD, idx[, .SD, .SDcols = -"idx"])]
}



#' Sample spatial points given a minimum distance between them
#'
#' @param dt [`icesat2.atl08_dt-class`] data.table. Which can be extracted using
#' [`ATL08_seg_attributes_dt()`]
#' @param size the sample size
#' @param radius the minimum radius between the points
#' @param spatialIndex the spatial index for accelerating the search
#' between points. default NULL will calculate a new ANN index tree.
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
sample <- function(dt, size, radius, spatialIndex, ...) {
  allMethods <- gsub("sample\\.", "", as.character(.S3methods("sample")))
  if (any(class(dt) %in% allMethods)) {
    UseMethod("sample")
  }

  base::sample(dt, size, ...)
}


saneSample <- function(x, n) {
  size_x <- length(x)
  if (n < 1 && n > 0) {
    return(sample(x, round(n * size_x)))
  }
  if (size_x < n) {
    warning(sprintf(
      "Not enough samples within group, truncated to %d samples",
      size_x
    ))
    n <- size_x
  }
  sample(x, n)
}

#' @export
`sample.icesat2.atl08_dt` <- radiusSamplingWorker


#' Sample spatial points given a minimum distance between them
#'
#' @param dt [`icesat2.atl03_seg_dt-class`] data.table. Which can be extracted using
#' [`ATL08_seg_attributes_dt()`]
#' @param size the sample size. Either an integer of absolute number of samples or
#' if it is between (0, 1) it will sample a percentage relative to the number of
#' available observations within the group.
#' @param radius the minimum radius between the points
#' @param spatialIndex the spatial index for accelerating the search
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
`sample.icesat2.atl03_seg_dt` <- radiusSamplingWorker




########################
## SAMPLING METHODS
########################
setRefClass("icesat2_sampling_method")

genericSamplingMethod <- function(fn, ...) {
  params <- list(...)
  params_string <- paste(gsub("^", ", ", names(params)), collapse = "")

  eval(parse(text = sprintf('
  function(size%1$s, ...) {
    methodContainer <- list(
      fn = fn,
      params = list(size%1$s, ...)
    )

    prepend_class(methodContainer, "icesat2_sampling_method")
    methodContainer
  }
  ', params_string)))
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
gridSampling <- genericSamplingMethod(gridSamplingWorker, grid_size = TRUE)

#' Get samples stratified by a variable binning histogram
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
stratifiedSampling <- genericSamplingMethod(stratifiedSamplingWorker, variable = TRUE)

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
radiusSampling <- genericSamplingMethod(radiusSamplingWorker, radius = TRUE)




# ggplot() +
#   geom_spatvector(data = geom, color = "orange", fill = NA) +
#   geom_point(data = sampled, aes(longitude, latitude, color = h_canopy), cex = 3) +
#   scale_x_continuous(breaks = seq(min(all_data$longitude), max(all_data$longitude), 1)) +
#   scale_y_continuous(breaks = seq(min(all_data$latitude), max(all_data$latitude), 1)) +
#   scale_color_gradient(low = "gray", high = "forestgreen") +
#   ggplot2::theme_bw()
