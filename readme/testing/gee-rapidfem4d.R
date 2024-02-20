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
    !(image[["Fmask"]] & cloudMask)
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

hls <- hls[["B2", "B3", "B4", "B5", "B6", "B7"]]
names(hls) <- c("blue", "green", "red", "nir", "swir1", "swir2")

# ndvi
hls[["ndvi"]] <- (hls[["nir"]] - hls[["red"]]) / (hls[["nir"]] + hls[["red"]])


# kndvi
sigma <- 1
nir <- hls[["nir"]]
red <- hls[["red"]]
nir_red <- nir - red
knr <- exp((nir_red)^2 / (2 * sigma^2))
hls[["kndvi"]] <- (1 - knr) / (1 + knr)
rng <- range(hls[["kndvi"]])

# evi
blue <- hls[["blue"]]
hls[["evi"]] <- 2.5 * (nir_red) / (nir + 6 * red - 7.5 * blue + 1)


# savi
hls[["savi"]] <- 1.5 * (nir_red) / (nir + red + 0.5)

# MSAVI2
p1 <- 2 * nir + 1
hls[["msavi"]] <- p1 - sqrt((p1^2) - 8 * (nir_red)) / 2
range(hls[["msavi"]])


# Linear spectral unmixing
soil <- c(0.14, 0.16, 0.22, 0.39, 0.45, 0.27)
veg <- c(0.086, 0.062, 0.043, 0.247, 0.109, 0.039)
water <- c(0.07, 0.039, 0.023, 0.031, 0.011, 0.007)

unmixing <- hls[[c("blue", "green", "red", "nir", "swir1", "swir2")]]
unmixing <- unmixing$unmix(list(soil, veg, water))
names(unmixing) <- c("f_soil", "f_veg", "f_water")
range(unmixing)
hls[["f_soil"]] <- unmixing[["f_soil"]]
hls[["f_veg"]] <- unmixing[["f_veg"]]
hls[["f_water"]] <- unmixing[["f_water"]]


# Simple Ratio Index
hls[["sri"]] <- nir / red

# NDWI
green <- hls[["green"]]
hls[["ndwi"]] <- (green - nir) / (green + nir)

# Green Chloropyl Index
hls[["gci"]] <- (nir / green) - 1

# WDRVI
hls[["wdrvi"]] <- ((0.1 * nir) - red) / ((0.1 * nir) + red)

# Global Vegetation Moisture Index
swir1 <- hls[["swir1"]]
hls[["gvmi"]] <- ((nir + 0.1) - (swir1 + 0.02)) / ((nir + 0.1) + (swir1 + 0.02))

# Chlorophyll Vegetation Index
hls[['cvi']] <- nir * (red / (green ^ 2))

# Clay minerals ratio
swir2 <- hls[["swir2"]]
hls[['cmr']] <- swir1 / swir2


# Radar vegetation index
hls[["rvi"]] <- sqrt(vv/(vv+vh))*(vv/vh)

img <- as.integer((hls[[c("blue", "green", "red", "nir", "swir1", "swir2")]]) * 1e4)
glcm <- img$glcmTexture(size = 3)

img <- hls$visualize(
  bands = c("cvi"),
  min = min(hls[["cvi"]]),
  max = max(hls[["cvi"]])
)



url <- img$getMapId()$tile_fetcher$url_format



library(leaflet)
library(magrittr)
# color_ramp <- leaflet::colorNumeric(palette = c("#FFFFFF", "#006400"), domain = c(0, 30))


leaflet_map <- leaflet::leaflet() %>%
  # addProviderTiles(providers$Esri.WorldImagery, group = "Other") %>%
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
