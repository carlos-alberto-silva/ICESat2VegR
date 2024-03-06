#' @title GridIndex Class
#'
#' @description
#' This class uses Grid Indexing system
#' for finding points that are within a specified range.
#'
#' @seealso
#' Mount, D. M.; Arya, S. ANN: A Library for Approximate Nearest Neighbor Searching,
#' available in <https://www.cs.umd.edu/~mount/ANN/>
#'
#' @export
GridIndex <- R6::R6Class("GridIndex",
  public = list(
    #' @field tree The C++ pointer for the built tree
    tree = NULL,
    #' @description
    #' Creates a new instance of the [`ANNIndex-class`] class.
    #'
    #' @param x NumericVector of x/longitude
    #' @param y NumericVector of y/latitude
    #' @param grid_size the size used for gridding
    initialize = function(x, y, grid_size) {
      self$tree <- pkg_module$GridIndex$new(x, y, grid_size)
    },
    #' @description
    #' Given a x, y point get the indexes that are within
    #' a specified radius.
    #'
    #' @param x NumericVector of x/longitude
    #' @param y NumericVector of y/latitude
    #' @param radius the search radius
    searchFixedRadius = function(x, y, radius) {
      self$tree$searchFixedRadius(x, y, radius)
    }
  ),
  cloneable = FALSE
)
