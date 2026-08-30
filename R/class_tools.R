#' Prepend a class to an object's list of classes
#'
#' @param obj The object to which prepend the class.
#' @param className [`character-class`] with the name of the class to prepend.
#' 
#' @return Nothing, it replaces the class attribute in place
#' 
#' @export
prepend_class <- function(obj, className) {
  call <- as.character(substitute(obj))
  parent <- parent.frame()
  obj2 <- parent[[call]]
  if (inherits(obj2, className)) {
    return()
  }

  attr(parent[[call]], "class") <- c(className, attr(parent[[call]], "class"))
}

#' Stop with an informative message if 'hdf5r' is not installed
#'
#' `hdf5r` is an optional (Suggests) dependency used only for reading and
#' writing local HDF5 (.h5) files.
#'
#' @return Nothing, called for the side effect of stopping if missing
#' @keywords internal
check_hdf5r <- function() {
  if (!requireNamespace("hdf5r", quietly = TRUE)) {
    stop(
      "Package 'hdf5r' is required to read/write local HDF5 (.h5) files but is not installed.\n",
      "Please install it with: install.packages(\"hdf5r\")",
      call. = FALSE
    )
  }
}
