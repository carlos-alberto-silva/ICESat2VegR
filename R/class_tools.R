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
