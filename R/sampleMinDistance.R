#' Sample spatial points given a minimum distance between them
#'
#' @param x NumericVector of x/longitude
#' @param y NumericVector of y/latitude
#' @param radius the minimum radius between the points
#' @param spatialIndex the spatial index for accelerating the search
#' between points. default NULL will calculate a new ANN index tree.
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
sampleMinDistance <- function(x, y, radius, nSamples, spatialIndex = NULL) {
  if (is.null(spatialIndex)) {
    spatialIndex <- ANNIndex$new(x, y)
  }

  return(pkg_module$sampleMinDistanceRcpp(x, y, radius, nSamples, spatialIndex$tree, as.integer(runif(1, 0, 2147483647))))
}
