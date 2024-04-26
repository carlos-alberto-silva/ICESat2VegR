#' @export
"[[<-.ee.ee_list.List" <- function(x, i, j, value) {
  x$add(value)
}


#' @export
`mean.ee.ee_list.List` <- function(x, ...) {
  if (!is.null(ee)) {
    return (x$reduce(ee$Reducer$mean()))
  }
}

#' @export
`max.ee.ee_list.List` <- function(..., na.rm = FALSE) {
  args <- list(...)
  x <- args[[1]]
  stopifnot("Expected only one argument for max for ee.List" = length(args) == 1)
  if (!is.null(ee)) {
    return (x$reduce(ee$Reducer$max()))
  }
}


#' @export
`min.ee.ee_list.List` <- function(x, na.rm = FALSE) {
  if (!is.null(ee)) {
    return (x$reduce(ee$Reducer$min()))
  }
}

