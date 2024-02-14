dt <- search_datasets("landsat 9", "tier 2", "toa")
(collection_id <- get_catalog_path(dt))

library(magrittr)
collection <- ee$ImageCollection(collection_id)

library(terra)
ext <- terra::ext(-87.6, -79.8, 24.5, 31)


filterBounds <- function(collection, ext) {
  bounds <- ee$Geometry$BBox(
    west = ext$xmin,
    south = ext$ymin,
    east = ext$xmax,
    north = ext$ymax
  )
  invisible(collection$filterBounds(bounds))
}

filterDate <- function(collection, ...) {
  args <- list(...)
  begin <- ifelse(length(args) == 2, args[[1]], args[[1]][1])
  end <- ifelse(length(args) == 2, args[[2]], args[[1]][2])

  invisible(collection$filterDate(begin, end))
}

setRefClass("ee.imagecollection.ImageCollection")
setMethod(
  "length",
  signature = c("ee.imagecollection.ImageCollection"),
  function(x) {
    return(x$size()$getInfo())
  }
)


collection %>%
  filterBounds(ext) %>%
  filterDate("2022-01-01", "2022-03-31") %>%
  length()


geometry$size()$getInfo()
extract_gee() <- function(points, img) {

}
