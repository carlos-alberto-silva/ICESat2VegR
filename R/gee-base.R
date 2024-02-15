setRefClass("ee.imagecollection.ImageCollection")
setRefClass("ee.image.Image")
#' Wraps around Google Earth Engine's filterDate method
#'
#' @param collection ee collection opened through
#' @param ... pass either a numeric vector with two ISO string dates
#' or two arguments with begin and end date range
#'
#' @export
filterDate <- function(collection, ...) {
  args <- list(...)
  begin <- ifelse(length(args) == 2, args[[1]], args[[1]][1])
  end <- ifelse(length(args) == 2, args[[2]], args[[1]][2])

  invisible(collection$filterDate(begin, end))
}

#' @export
setMethod(
  "length",
  signature = c("ee.imagecollection.ImageCollection"),
  function(x) {
    return(x$size()$getInfo())
  }
)

#' Get band names from a google earth engine collection
#' @export
bandNames <- function(x) {
  x$getMapId()$image$bandNames()$getInfo()
}

#' @export
setMethod(
  "median",
  signature = c("ee.imagecollection.ImageCollection"),
  function(x) {
    invisible(x$median())
  }
)

#' @export
getTileUrl <- function(img) {
  img$getMapId()$tile_fetcher$url_format
}

#' @export
filter <- dplyr::filter

#' @export
setMethod(
  "filter",
  signature = c("ee.imagecollection.ImageCollection"),
  function(.data, ...) {
    args <- list(...)
    invisible(.data$filter(args[[1]]))
  }
)

#' @export
map <- function(collection, fn) {
  useMethod("map")
}

#' @export
setMethod(
  "map",
  signature = c("ee.imagecollection.ImageCollection"),
  function(collection, fn) {
    collection$map(fn)
  }
)



#' @export
"*.ee.image.Image" <- function(e1, e2) e1$multiply(e2)

#' @export
"/.ee.image.Image" <- function(e1, e2) e1$divide(e2)

#' @export
"+.ee.image.Image" <- function(e1, e2) e1$add(e2)

#' @export
"-.ee.image.Image" <- function(e1, e2) e1$subtract(e2)

#' @export
"|.ee.image.Image" <- function(e1, e2) e1$bitwiseOr(e2)

#' @export
"&.ee.image.Image" <- function(e1, e2) e1$bitwiseAnd(e2)

#' @export
"!.ee.image.Image" <- function(e1) e1$bitwiseNot()

#' @export
"^.ee.image.Image" <- function(e1, e2) e1$pow(e2)

#' @export
"exp.ee.image.Image" <- function(e1) e1$exp()

#' @export
"log.ee.image.Image" <- function(e1) e1$log()

#' @export
">.ee.image.Image" <- function(e1, e2) e1$gt(e2)

#' @export
"<.ee.image.Image" <- function(e1, e2) e1$lt(e2)

#' @export
">=.ee.image.Image" <- function(e1, e2) e1$gte(e2)

#' @export
"<=.ee.image.Image" <- function(e1, e2) e1$lte(e2)

#' @export
"%%.ee.image.Image" <- function(e1, e2) e1$mod(e2)

#' @export
"as.integer.ee.image.Image" <- function(e1) e1$int()

#' @export
"[[.ee.image.Image" <- function(x, ...) do.call(x$select, list(...))

#' @export 
"names.ee.image.Image" <- function(x) x$bandNames()$getInfo()

#' @export 
"names<-.ee.image.Image" <- function(x, value) x$select(names(x), value)

#' @export 
"[[<-.ee.image.Image" <- function(x, i, j, ..., value) {
  value2 <- value$rename(i)
  x$addBands(value2)
}
