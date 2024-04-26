#' Prepend a class to an object's list of classes
#'
#' @param obj The object to which prepend the class.
#' @param className [`character-class`] with the name of the class to prepend.
#' 
#' @return Nothing, it replaces the class attribute in place
#' 
#' @import lazyeval
#' @export
prepend_class <- function(obj, className) {
  call <- lazyeval::expr_text(obj)
  parent <- parent.frame()
  obj2 <- parent[[call]]
  if (inherits(obj2, className)) {
    return()
  }

  attr(parent[[call]], "class") <- c(className, attr(parent[[call]], "class"))
}
