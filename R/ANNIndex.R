#' @title ANNIndex Class
#'
#' @description
#' This class uses Approximate Neareast Neighbor Index
#' for finding points that are within a specified range.
#'
#' @seealso
#' Mount, D. M.; Arya, S. ANN: A Library for Approximate Nearest Neighbor Searching,
#' available in <https://www.cs.umd.edu/~mount/ANN/>
#'
#' @export
ANNIndex <- R6::R6Class("ANNIndex",
  public = list(
    #' @field tree The C++ pointer for the built tree
    tree = NULL,
    #' @description
    #' Creates a new instance of the [`ANNIndex`] class.
    #'
    #' @param x NumericVector of x/longitude
    #' @param y NumericVector of y/latitude
    #' @param radius the minimum radius between the points
    initialize = function(x, y) {
      self$tree <- pkg_module$ANNIndex$new(x, y)
    },
    #' @description
    #' Given a x, y point get the indexes that are within
    #' a specified radius.
    #'
    #' @param x NumericVector of x/longitude
    #' @param y NumericVector of y/latitude
    #' @param radius the minimum radius between the points
    searchFixedRadius = function(x, y, radius) {
      self$tree$searchFixedRadius(x, y, radius)
    }
  ),
  cloneable = FALSE
)
