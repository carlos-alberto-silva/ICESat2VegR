#' @export
"[[.ee.feature.Feature" <- function(x, i, j) x$getNumber(i)


#' @export
"[[<-.ee.feature.Feature" <- function(x, i, j, value) {
  the_list <- list()
  the_list[[i]] <- value
  x$set(the_list)
}
