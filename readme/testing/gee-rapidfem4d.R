source("R/gee-search.R")
source("R/gee-base.R")
geemap <- reticulate::import("geemap")

library(ICESat2VegR)
library(terra)
library(magrittr)

ext <- terra::ext(-87.6, -79.8, 24.5, 31)
aoi <- ee$Geometry$BBox(
  west = ext$xmin,
  south = ext$ymin,
  east = ext$xmax,
  north = ext$ymax
)

search <- search_datasets("hls", "landsat")

# get_catalog_id <- get_catalog_path
id <- get_catalog_id(search)
collection <- ee$ImageCollection(id)


cloudMask <-
  2^0 +
  2^1 +
  2^2

hlsMask <- function(image) {
  image$updateMask(
    image[["Fmask"]] & cloudMask
  )
}

waterMask <- function(image) {
  image$updateMask(
    image[["B5"]] >= 0.2
  )
}


hls <- collection$
  filterBounds(aoi)$
  filterDate("2021-02-01", "2021-05-31")$
  filter("CLOUD_COVERAGE < 10")$
  map(hlsMask)$
  map(waterMask)$
  median()

hls <- as.integer(hls * 10000)
hls <- hls[['B2', 'B3', 'B4', 'B5', 'B6','B7']]
names(hls) <- c('blue', 'green','red', 'nir', 'swir1', 'swir2')
hls[['ndvi']] <- (hls[['nir']] - hls[['red']]) / (hls[['nir']] + hls[['red']])

 D2 = im.select('nir').subtract(im.select('red')).pow(2).rename('d2')
    gamma = ee.Number(4e6).multiply(-2.0)
    k = D2.divide(gamma).exp()
    kndvi = ee.Image.constant(1).subtract(k).divide(ee.Image.constant(1).add(k)).rename('kndvi')
    return im.addBands(kndvi).multiply(1000).int()


img <- hls$visualize(
  bands = c("ndvi"),
  min = -1,
  max = 1
)


url <- img$getMapId()$tile_fetcher$url_format



library(leaflet)
library(magrittr)
# color_ramp <- leaflet::colorNumeric(palette = c("#FFFFFF", "#006400"), domain = c(0, 30))


leaflet_map <- leaflet::leaflet() %>%
addProviderTiles(providers$Esri.WorldImagery, group = "Other") %>%
  leaflet::addTiles(
    urlTemplate = url,
    options = leaflet::tileOptions(opacity = 1),
    group = "Landsat"
  ) %>%
  addLayersControl(
    overlayGroups = c("Landsat"),
    options = layersControlOptions(collapsed = FALSE)
  ) %>%
  leaflet::setView(lng = -82, lat = 28, zoom = 11)

leaflet_map
