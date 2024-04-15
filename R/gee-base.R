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
    return(invisible(x$size()$getInfo()))
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
    return(invisible(collection$map(fn)))
  }
)



#' @export
"*.ee.image.Image" <- function(e1, e2) {
  if (!inherits(e1, "ee.image.Image")) {
    return(invisible(ee$Image$constant(e1)$multiply(e2)))
  }
  return(invisible(e1$multiply(e2)))
}

#' @export
"/.ee.image.Image" <- function(e1, e2) {
  if (!inherits(e1, "ee.image.Image")) {
    return(invisible(ee$Image$constant(e1)$divide(e2)))
  }
  return(invisible(e1$divide(e2)))
}

#' @export
"+.ee.image.Image" <- function(e1, e2) {
  if (!inherits(e1, "ee.image.Image")) {
    return(invisible(ee$Image$constant(e1)$add(e2)))
  }
  return(invisible(e1$add(e2)))
}

#' @export
"-.ee.image.Image" <- function(e1, e2) {
  if (!inherits(e1, "ee.image.Image")) {
    return(invisible(ee$Image$constant(e1)$subtract(e2)))
  }
  return(invisible(e1$subtract(e2)))
}

#' @export
"|.ee.image.Image" <- function(e1, e2) {
  if (!inherits(e1, "ee.image.Image")) {
    return(invisible(ee$Image$constant(e1)$bitwiseOr(e2)))
  }
  return(invisible(e1$bitwiseOr(e2)))
}

#' @export
"&.ee.image.Image" <- function(e1, e2) {
  if (!inherits(e1, "ee.image.Image")) {
    return(invisible(ee$Image$constant(e1)$bitwiseAnd(e2)))
  }
  return(invisible(e1$bitwiseAnd(e2)))
}

#' @export
"!.ee.image.Image" <- function(e1) e1$Not()

#' @export
"^.ee.image.Image" <- function(e1, e2) {
  if (!inherits(e1, "ee.image.Image")) {
    return(invisible(ee$Image$constant(e1)$pow(e2)))
  }
  return(invisible(e1$pow(e2)))
}

#' @export
"exp.ee.image.Image" <- function(e1) {
  return(invisible(e1$exp()))
}

#' @export
"log.ee.image.Image" <- function(e1) {
  return(invisible(e1$log()))
}

#' @export
">.ee.image.Image" <- function(e1, e2) {
  if (!inherits(e1, "ee.image.Image")) {
    return(invisible(ee$Image$constant(e1)$gt(e2)))
  }
  return(invisible(e1$gt(e2)))
}

#' @export
"<.ee.image.Image" <- function(e1, e2) {
  if (!inherits(e1, "ee.image.Image")) {
    return(invisible(ee$Image$constant(e1)$lt(e2)))
  }
  return(invisible(e1$lt(e2)))
}

#' @export
">=.ee.image.Image" <- function(e1, e2) {
  if (!inherits(e1, "ee.image.Image")) {
    return(invisible(ee$Image$constant(e1)$gte(e2)))
  }
  return(invisible(e1$gte(e2)))
}

#' @export
"<=.ee.image.Image" <- function(e1, e2) {
  if (!inherits(e1, "ee.image.Image")) {
    return(invisible(ee$Image$constant(e1)$lte(e2)))
  }
  return(invisible(e1$lte(e2)))
}

#' @export
"%%.ee.image.Image" <- function(e1, e2) {
  if (!inherits(e1, "ee.image.Image")) {
    return(invisible(ee$Image$constant(e1)$mod(e2)))
  }
  return(invisible(e1$mod(e2)))
}

#' @export
"as.integer.ee.image.Image" <- function(e1) {
  return(invisible(e1$int()))
}

#' @export
"[[.ee.image.Image" <- function(x, ...) {
  return(invisible(do.call(x$select, list(...))))
}

#' @export
"names.ee.imagecollection.ImageCollection" <- function(x) {
  return(invisible(x$first()$bandNames()$getInfo()))
}


#' @export
"names.ee.image.Image" <- function(x) {
  return(invisible(x$bandNames()$getInfo()))
}

#' @export
"names<-.ee.image.Image" <- function(x, value) {
  return(invisible(x$select(names(x), value)))
}

#' @export
"[[<-.ee.image.Image" <- function(x, i, j, ..., value) {
  if (is.character(i) && i != "") {
    value2 <- value$rename(i)
  } else {
    value2 <- value
  }
  return(invisible(x$addBands(value2)))
}

#' @export
"mean.ee.image.Image" <- function(x) {
  return(invisible(x$getMapId()$tile_fetcher$url_format))
}

#' @export
"min.ee.image.Image" <- function(x, ...) {
  args <- list(...)
  args[["na.rm"]] <- NULL
  args$reducer <- ee$Reducer$min()
  if (is.null(args$geometry)) {
    args$geometry <- ee$Geometry$BBox(
      west = -180,
      south = -90,
      east = 180,
      north = 90
    )
  }
  return(invisible(unlist(
    do.call(x$reduceRegion, args)$
      getInfo()
  )))
}

#' @export
"max.ee.image.Image" <- function(x, ...) {
  args <- list(...)
  args[["na.rm"]] <- NULL
  args$reducer <- ee$Reducer$max()
  if (is.null(args$geometry)) {
    args$geometry <- ee$Geometry$BBox(
      west = -180,
      south = -90,
      east = 180,
      north = 90
    )
  }
  return(invisible(unlist(
    do.call(x$reduceRegion, args)$
      getInfo()
  )))
}

#' @export
"range.ee.image.Image" <- function(x, ...) {
  return(invisible(rev(unlist(x$reduceRegion(
    reducer = ee$Reducer$minMax(),
    geometry = ee$Geometry$BBox(
      west = -180,
      south = -90,
      east = 180,
      north = 90
    )
  )$getInfo()))))
}

#' @export
"sqrt.ee.image.Image" <- function(x, ...) {
  return(invisible(x$sqrt()))
}


#' @export
glcmTexture <- function(x, ...) {
  UseMethod("glcmTexture")
}

#' @export
"glcmTexture.ee.image.Image" <- function(x, ...) {
  return(invisible(x$glcmTexture(...)))
}

#' @export
unmix <- function(x, ...) {
  UseMethod("unmix")
}


#' @export
"unmix.ee.image.Image" <- function(x, ...) {
  return(invisible(x$unmix(...)))
}


#' @export
"c.ee.imagecollection.ImageCollection" <- function(x, ...) {
  args <- list(...)
  result <- x
  for (arg in args) {
    result <- result$merge(arg)
  }
  return(invisible(result))
}


#' @export
"c.ee.image.Image" <- function(x, ...) {
  return(invisible(ee$Image$cat(x, ...)))
}


#' @export
"print.ee.image.Image" <- function(x, ...) {
  cat("ee.image.Image\n\n")
  cat("Bands\n")
  print(names(x))
  invisible()
}


#' @export
slope <- function(x, ...) {
  UseMethod("slope")
}

#' @export
"slope.ee.image.Image" <- function(x, ...) {
  return(invisible(ee$Terrain$slope(x, ...)))
}

#' @export
aspect <- function(x, ...) {
  UseMethod("aspect")
}

#' @export
"aspect.ee.image.Image" <- function(x, ...) {
  return(invisible(ee$Terrain$aspect(x, ...)))
}
