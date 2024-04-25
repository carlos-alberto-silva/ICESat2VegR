setRefClass("ee.imagecollection.ImageCollection")
setRefClass("ee.image.Image")


#' @export
setMethod(
  "length",
  signature = c("ee.imagecollection.ImageCollection"),
  function(x) {
    return(invisible(x$size()$getInfo()))
  }
)

#' Retrieve Google Earth Engine's tile url for an Image or ImageCollection
#'
#' @param img The `Image` or `ImageCollection` to retrieve the URL
#'
#' @return The url for the tile service.
#'
#' @export
getTileUrl <- function(img) {
  img$getMapId()$tile_fetcher$url_format
}

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
  return(invisible(x$addBands(value2, overwrite = TRUE)))
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


#' Maps to `ee.Image.glcmTexture`
#'
#' Computes texture metrics from the Gray Level Co-occurrence Matrix around
#' each pixel of every band. The GLCM is a tabulation of how often different
#' combinations of pixel brightness values (grey levels) occur in an image. It
#' counts the number of times a pixel of value X lies next to a pixel of value Y, in a
#' particular direction and distance. and then derives statistics from this
#' tabulation.
#'
#' @param x The input image to calculate the texture on.
#' @param size integer, default 1. The size of the neighborhood to include in each GLCM.
#' @param kernel default NULL. A kernel specifying the x and y offsets over which to compute
#' the GLCMs. A GLCM is computed for each pixel in the kernel that is non-zero, except the
#' center pixel and as long as a GLCM hasn't already been computed for the same direction and
#' distance. For example, if either or both of the east and west pixels are set, only 1
#' (horizontal) GLCM is computed. Kernels are scanned from left to right and top to bottom.
#' The default is a 3x3 square, resulting in 4 GLCMs with the offsets
#' (-1, -1), (0, -1), (1, -1) and (-1, 0).
#' @param average logical, default TRUE. If true, the directional bands for each metric are
#' averaged.


#'
#' @return Another Earth Engine image with bands described on details section.
#'
#' @details
#'
#' This implementation computes the 14 GLCM metrics proposed by Haralick, and
#' 4 additional metrics from Conners. Inputs are required to be integer valued.
#'
#' The output consists of 18 bands per input band if directional averaging is on
#' and 18 bands per directional pair in the kernel, if not:
#'
#' 1. ASM: f1, Angular Second Moment; measures the number of repeated pairs
#' 1. CONTRAST: f2, Contrast; measures the local contrast of an image
#' 1. CORR: f3, Correlation; measures the correlation between pairs of pixels
#' 1. VAR: f4, Variance; measures how spread out the distribution of gray-levels is
#' 1. IDM: f5, Inverse Difference Moment; measures the homogeneity
#' 1. SAVG: f6, Sum Average
#' 1. SVAR: f7, Sum Variance
#' 1. SENT: f8, Sum Entropy
#' 1. ENT: f9, Entropy. Measures the randomness of a gray-level distribution
#' 1. DVAR: f10, Difference variance
#' 1. DENT: f11, Difference entropy
#' 1. IMCORR1: f12, Information Measure of Corr. 1
#' 1. IMCORR2: f13, Information Measure of Corr. 2
#' 1. MAXCORR: f14, Max Corr. Coefficient. (not computed)
#' 1. DISS: Dissimilarity
#' 1. INERTIA: Inertia
#' 1. SHADE: Cluster Shade
#' 1. PROM: Cluster prominence
#'
#' @seealso
#' More information can be found in the two papers: Haralick et. al, 'Textural
#' Features for Image Classification', http://doi.org/10.1109/TSMC.1973.4309314
#' and Conners, et al, Segmentation of a high-resolution urban scene using texture
#' operators', http://doi.org/10.1016/0734-189X(84)90197-X.
#'
#' https://developers.google.com/earth-engine/apidocs/ee-image-glcmtexture
#'
#' @export
glcmTexture <- function(x, size = 1, kernel = NULL, average = TRUE) {
  UseMethod("glcmTexture")
}

#' @export
"glcmTexture.ee.image.Image" <- function(x, size = 1, kernel = NULL, average = TRUE) {
  return(invisible(x$glcmTexture(size = size, kernel = kernel, average = average)))
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

#' Calculates slope in degrees from a terrain DEM.
#'
#' The local gradient is computed using the 4-connected neighbors of each pixel,
#' so missing values will occur around the edges of an image.
#'
#' @param x The `ee.Image` on which to apply the aspect function.
#'
#' @return An `ee.Image` with a single band named "Slope".
#'
#' @export
slope <- function(x) {
  UseMethod("slope")
}


#' @export
"slope.ee.image.Image" <- function(x) {
  return(invisible(ee$Terrain$slope(x)))
}

#' Calculates aspect in degrees from a terrain DEM.
#'
#' The local gradient is computed using the 4-connected neighbors of each
#' pixel, so missing values will occur around the edges of an image.
#'
#' @param x The `ee.Image` on which to apply the aspect function.
#'
#' @return An `ee.Image` with a single band named "Aspect".
#'
#' @export
aspect <- function(x) {
  UseMethod("aspect")
}

#' @export
"aspect.ee.image.Image" <- function(x) {
  return(invisible(ee$Terrain$aspect(x)))
}
