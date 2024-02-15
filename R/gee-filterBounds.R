
setRefClass("ee.imagecollection.ImageCollection")


#' Wraps around Google Earth Engine's filterBound method
#' @export
filterBounds <- function(
    collection,
    bounds,
    ...) {
  UseMethod("filterBounds")
}

#' @importClassesFrom terra SpatExtent
setMethod(
  "filterBounds",
  signature = c(
    "ee.imagecollection.ImageCollection",
    "SpatExtent"
  ),
  function(collection, bounds) {
    bounds <- ee$Geometry$BBox(
      west = bounds$xmin,
      south = bounds$ymin,
      east = bounds$xmax,
      north = bounds$ymax
    )
    invisible(collection$filterBounds(bounds))
  }
)

#' @importClassesFrom terra SpatVector
setMethod(
  "filterBounds",
  signature = c(
    "ee.imagecollection.ImageCollection",
    "SpatVector"
  ),
  function(collection, bounds) {
    ext <- terra::ext(bounds)
    bounds <- ee$Geometry$BBox(
      west = ext$xmin,
      south = ext$ymin,
      east = ext$xmax,
      north = ext$ymax
    )
    invisible(collection$filterBounds(bounds))
  }
)

setMethod(
  "filterBounds",
  signature = c(
    "ee.imagecollection.ImageCollection",
    "numeric"
  ),
  function(collection, bounds, ...) {
    args <- list(...)
    stopifnot("expected a numeric vector or four arguments for binding" = {
      length(args) == 3
    })
    bounds <- ee$Geometry$BBox(
      west = bounds$xmin,
      south = bounds$ymin,
      east = bounds$xmax,
      north = bounds$ymax
    )
    invisible(collection$filterBounds(bounds))
  }
)
