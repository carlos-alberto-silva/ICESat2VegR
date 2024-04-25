source("R/gee-search.R")
source("R/gee-class.R")

# require(pacman)
# p_load(reticulate)
# conda_create(envname = "ICEsat/2VegR", path = "C:/Users/c.silva/AppData/Local/r-miniconda/envs/")
# use_condaenv("C:/Users/c.silva/AppData/Local/r-miniconda/envs/ICEsat-2VegR", required = TRUE)

# py_run_string("import sys; print(sys.version)")

# Sys.setenv(
#   PATH = paste(Sys.getenv("PATH"),
#     "C:/Users/c.silva/AppData/Local/r-miniconda/envs/ICEsat-2VegR/DLLs",
#     sep = ";"
#   )
# )

# reticulate::install_miniconda()
# reticulate::py_install(c("earthengine-api", "earthaccess", "h5py"))


ee <- reticulate::import("ee")
geemap <- reticulate::import("geemap")

geemap$ee_initialize()
library(dplyr)

dt <- search_datasets('landsat 9', 'tier 2','toa')
(collection_id <- get_catalog_path(dt))

collection <- 
  ee$ImageCollection(collection_id)


bounds <- ee$Geometry$BBox(
  west = -87.6,
  south = 24.5,
  east = -79.8,
  north = 31
  )

collection2 <- collection$filterBounds(bounds)$filterDate("2022-01-01", "2022-01-31")
collection2$size()$getInfo()
mapId <- collection2$getMapId()
url <- mapId$tile_fetcher$url_format

library(leaflet)
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


map = geemap$Map()
map$search_loc_geom



feature_list[[1]]$image$bandNames()$getInfo()
expression  = "landsat.B2"

roi <- ee$Geometry$Rectangle(c(-87.6, 24.5, -79.8, 31))
l8 <- ee$ImageCollection("LANDSAT/LC08/C01/T1")





image <- ee$Algorithms$Landsat$simpleComposite(
  collection = l8$filterDate("2018-01-01", "2018-01-31"), asFloat = TRUE
)




library(ICESat2VegR)
atl08_h5 <- ICESat2VegR::ATL08_read("../inst/exdata/ATL08_20190819062101_07970402_006_02.h5")
atl08_seg_dt <- ATL08_seg_attributes_dt(
  atl08_h5,
  beam = c("gt1l", "gt2l", "gt3l"),
  attribute = c(
    "h_canopy",
    "h_te_mean",
    "night_flag",
    "segment_landcover"
  )
)
dataset <- atl08_seg_dt[night_flag == 1 & !segment_landcover %in% c(0, 80, 200) & h_canopy < 50]
set.seed(5789174)
shuffled <- dataset[order(runif(nrow(dataset)))]

train <- shuffled[1:700]
test <- shuffled[701:1000]

library(terra)
v <- terra::vect(
  as.data.frame(train[, list(latitude, longitude, h_canopy)]),
  geom = c("longitude", "latitude"),
  crs = "epsg:4326"
)

terra::writeVector(v, "../train.geojson", filetype = "geojson", overwrite = TRUE)


v2 <- terra::vect(
  as.data.frame(test[, list(latitude, longitude, h_canopy)]),
  geom = c("longitude", "latitude"),
  crs = "epsg:4326"
)

terra::writeVector(v2, "../test.geojson", filetype = "geojson", overwrite = TRUE)


# geojson <- geemap$geojson_to_ee("../train.geojson")

# sampled <- image$sampleRegions(
#   collection = geojson,
#   scale = image$projection()$nominalScale()
# )

library(data.table)
backJson <- sampled$getInfo()
dt <- rbindlist(
  lapply(
    backJson$features,
    function(x) {
      dt <- as.data.table(x$properties)
      dt[, id := gsub("_\\d+", "", x$id)]
      dt
    }
  )
)
dt[, h_canopy := NULL]

v <- terra::vect(
  cbind(train, dt),
  geom = c("longitude", "latitude"),
  crs = "epsg:4326"
)

terra::writeVector(v, "../train.geojson", filetype = "geojson", overwrite = TRUE)

trained_classifier <- ee$Classifier$smileRandomForest(100)$setOutputMode("REGRESSION")
trained_classifier <- trained_classifier$train(
  features = sampled,
  classProperty = "h_canopy",
  inputProperties = feature_names,
)

# Define Florida Bounding Box
roi <- ee$Geometry$Rectangle(c(-87.6, 24.5, -79.8, 31))


classified2 <- image$select(feature_names)$clip(roi)$classify(trained_classifier)
# Define the visualization parameters


color_image <- classified2$visualize(
  bands = "classification",
  min = 0,
  max = 30,
  palette = c("FFFFFF", "006400")
)


map_id <- ee$data$getMapId(list(image = color_image))
url <- map_id[["tile_fetcher"]]$url_format

library(leaflet)
color_ramp <- leaflet::colorNumeric(palette = c("#FFFFFF", "#006400"), domain = c(0, 30))

leaflet_map <- leaflet::leaflet() %>%
  addProviderTiles(providers$Esri.WorldImagery, group = "Other") %>%
  leaflet::addTiles(
    urlTemplate = url,
    options = leaflet::tileOptions(opacity = 1),
    group = "Classification"
  ) %>%
  addLayersControl(
    overlayGroups = c("Classification"),
    options = layersControlOptions(collapsed = FALSE)
  ) %>%
  leaflet::setView(lng = -82, lat = 28, zoom = 11) %>%
  addLegend(
    title = "h_canopy",
    position = "bottomright",
    pal = color_ramp,
    values = c(0, 5, 10, 15, 20, 25, 30), # Adjust breaks as needed
    opacity = 1
  )

leaflet_map
