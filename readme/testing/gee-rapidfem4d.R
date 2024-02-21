# source("R/gee-search.R")
# source("R/gee-base.R")
# geemap <- reticulate::import("geemap")


library(ICESat2VegR)
library(terra)
library(magrittr)

geom <- terra::vect("../output.gpkg")
ext <- terra::ext(geom)
center <- terra::geom(terra::centroids(terra::vect(ext)))[, c("x", "y")]
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

# evi
blue <- hls[["blue"]]
hls[["evi"]] <- 2.5 * (nir_red) / (nir + 6 * red - 7.5 * blue + 1)


# savi
hls[["savi"]] <- 1.5 * (nir_red) / (nir + red + 0.5)

# MSAVI2
p1 <- 2 * nir + 1
hls[["msavi"]] <- p1 - sqrt((p1^2) - 8 * (nir_red)) / 2


# Linear spectral unmixing
soil <- c(0.14, 0.16, 0.22, 0.39, 0.45, 0.27)
veg <- c(0.086, 0.062, 0.043, 0.247, 0.109, 0.039)
water <- c(0.07, 0.039, 0.023, 0.031, 0.011, 0.007)

img <- hls[[c("blue", "green", "red", "nir", "swir1", "swir2")]]
unmixing <- img$unmix(list(soil, veg, water))

names(unmixing) <- c("f_soil", "f_veg", "f_water")


hls[[""]] <- unmixing
# Or
# combination <- c(hls, unmixing)

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
hls[["cvi"]] <- nir * (red / (green^2))

# Clay minerals ratio
swir2 <- hls[["swir2"]]
hls[["cmr"]] <- swir1 / swir2


# Texture
img <- as.integer((hls[[c("blue", "green", "red", "nir", "swir1", "swir2")]]) * 1e4)
hls[[""]] <- img$glcmTexture(size = 3)


# Radar vegetation index
# hls[["rvi"]] <- sqrt(vv/(vv+vh))*(vv/vh)

###########################
# Elevation
###########################
result <- search_datasets("dem", "copernicus")
result

catalog_id <- get_catalog_id(result)
catalog_id

elevation <- ee$ImageCollection(catalog_id)
elevation <- elevation$
  filterBounds(aoi)$
  median()
elevation <- elevation$clip(aoi)
elevation


the_slope <- slope(elevation)
the_aspect <- elevation %>% aspect()

# stackDem
stackDem <- c(elevation, the_slope, the_aspect)


stackS2 <- hls[["ndvi", "kndvi", "evi", "savi", "msavi", "sri", "ndwi", "gci", "wdrvi", "gvmi", "cvi", "cmr"]]
stackS2[[""]] <- img
stackS2[[""]] <- unmixing

fullStack <- c(hls, stackDem)

kernel <- ee$Kernel$fixed(3, 3, list(c(1, 1, 1), c(1, 1, 1), c(1, 1, 1)), 3, 3, FALSE)
reducerNames <- c("mean", "min", "max", "stdDev")
for (reducerName in reducerNames) {
  reducer <- ee$Reducer[[reducerName]]()
  fullStack[[""]] <- stackS2$reduceNeighborhood(reducer, kernel)
  fullStack[[""]] <- stackDem$reduceNeighborhood(reducer, kernel)
}
fullStack



visBand <- "ndvi_mean"

img <- fullStack$visualize(
  bands = c(visBand),
  min = min(fullStack[[visBand]]),
  max = max(fullStack[[visBand]])
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
  leaflet::setView(lng = center[1], lat = center[2], zoom = 11)

leaflet_map





# Extract
geom <- sf::st_read("../output.gpkg")
sample_size <- 1000


geom <- geom %>%
  sf::st_sample(size = sample_size) %>%
  sf::st_cast("POINT") %>%
  sf::st_sf()

if (file.exists("../output.geojson")) {
  unlink("../output.geojson")
}

geom %>% sf::st_write("../output.geojson")


parsed <- jsonlite::parse_json(readLines("../output.geojson"))

geojson <- ee$FeatureCollection(parsed)

sampled <- fullStack$sampleRegions(
  collection = geojson,
  scale = fullStack$projection()$nominalScale()
)


downloadId <- ee$data$getTableDownloadId(list(table = sampled, format = "csv"))
url2 = ee$data$makeTableDownloadUrl(downloadId)
res <- ee$data$requests$get(url2)
df <- read.table(
  text = res$content$decode('utf8'),
  sep = ',',
  header = TRUE
)
ids <- df[,1]
ids <- as.integer(gsub("_0", "", ids)) + 1
df[,1] <- ids

coords <- sf::st_coordinates(geom)[ids, ]

final_vect <- terra::vect(
  cbind(coords, df),
  geom = c("X", "Y"),
  crs = "epsg:4326"
)

terra::writeVector(final_vect, "../final_output.gpkg")
