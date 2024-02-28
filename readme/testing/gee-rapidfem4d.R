# source("R/gee-search.R")
# source("R/gee-base.R")
# geemap <- reticulate::import("geemap")


library(ICESat2VegR)
library(terra)
library(magrittr)

geom <- terra::vect("../train.geojson")
geom2 <- terra::vect("../test.geojson")
all_geoms <- terra::union(geom, geom2)
ext <- terra::ext(all_geoms)
center <- as.numeric(terra::geom(terra::centroids(terra::vect(ext)))[, c("x", "y")])
aoi <- ee$Geometry$BBox(
  west = ext$xmin,
  south = ext$ymin,
  east = ext$xmax,
  north = ext$ymax
)

aoi <- aoi$buffer(30)

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
  filterDate("2020-02-01", "2020-05-31")$
  filter("CLOUD_COVERAGE < 10")$
  map(hlsMask)$
  map(waterMask)$
  median()



hls <- hls[["B2", "B3", "B4", "B5", "B6", "B7"]]
names(hls) <- c("blue", "green", "red", "nir", "swir1", "swir2")
hls <- hls$clip(aoi)

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
result <- search_datasets("nasa", "dem")
catalog_id <- get_catalog_id(result)

elevation <- ee$Image(catalog_id)

the_slope <- as.integer(slope(as.integer(elevation)) * 1000)
the_aspect <- elevation %>% aspect()

# stackDem
stackDem <- c(elevation, the_slope, the_aspect)$clip(aoi)

fullStack <- c(hls, stackDem)

kernel <- ee$Kernel$fixed(3, 3, list(c(1, 1, 1), c(1, 1, 1), c(1, 1, 1)), 3, 3, FALSE)
reducerNames <- c("mean", "min", "max", "stdDev")
for (reducerName in reducerNames) {
  reducer <- ee$Reducer[[reducerName]]()
  fullStack[[""]] <- stackS2$reduceNeighborhood(reducer, kernel)
  fullStack[[""]] <- stackDem$reduceNeighborhood(reducer, kernel)
}
fullStack


source("R/gee-base.R")

visBand <- "slope"
img <- fullStack[[visBand]]

img <- img$visualize(
  bands = visBand,
  min = min(img, scale = 900, geometry = aoi),
  max = max(img, scale = 900, geometry = aoi)
)

url <- img$getMapId()$tile_fetcher$url_format
library(leaflet)
library(magrittr)
# color_ramp <- leaflet::colorNumeric(palette = c("#FFFFFF", "#006400"), domain = c(0, 30))
coords <- terra::geom(geom)

leaflet_map <- leaflet::leaflet() %>%
  addProviderTiles(providers$Esri.WorldImagery, group = "Other") %>%
  leaflet::addTiles(
    urlTemplate = url,
    options = leaflet::tileOptions(opacity = 1),
    group = "Landsat"
  ) %>%
  leaflet::addCircleMarkers(
    lng = coords[, "x"],
    lat = coords[, "y"],
    radius = 2,
    stroke = FALSE,
    fillOpacity = 1,
    fillColor = "yellow"
  ) %>%
  addLayersControl(
    overlayGroups = c("Landsat"),
    options = layersControlOptions(collapsed = FALSE)
  ) %>%
  leaflet::setView(lng = center[1], lat = center[2], zoom = 11)

leaflet_map





source("R/extract.R")

# Extract

train <- terra::vect("../train.geojson")
test <- terra::vect("../test.geojson")

train_extract <- extract(fullStack, train, 30)
test_extract <- extract(fullStack, test, 30)
train_dt <- ee_to_dt(train_extract)
test_dt <- ee_to_dt(test_extract)

x_columns <- names(train_dt)[-c(1, which(names(train_dt) == "h_canopy"))]
x <- train_dt[, .SD, .SDcols = x_columns]
y <- train_dt$h_canopy
x_test <- test_dt[, .SD, .SDcols = x_columns]
y_test <- test_dt$h_canopy

randomForestClassifier <- randomForestRegression(train_extract, property_name = "h_canopy", nTrees = 100, nodesize = 1)
gee_predict <- randomForestClassifier %>% predict(test_extract)

StatModel(y_test, gee_predict2, xlim = c(0, 30), ylim = c(0, 30))

(gee_rmse <- sqrt(mean((y_test - gee_predict)^2)))
lims <- c(0, max(c(y_test, gee_predict)))
plot(y_test, gee_predict, xlim = lims, ylim = lims, main = sprintf("GEE randomForest (RMSE = %.4f)", gee_rmse))
lines(c(-9999, 9999), c(-9999, 9999), col = "red")

library(randomForest)
rf <- randomForest::randomForest(x, y, ntree = 100, nodesize = 1)
rf_predict <- rf %>% predict(test_dt)
(rf_rmse <- sqrt(mean((y_test - rf_predict)^2)))
trees <- build_forest(rf)


gee_rf <- ee$Classifier$decisionTreeEnsemble(trees)
gee_rf_predict <- gee_rf %>% predict(test_extract)
(gee_rf_rmse <- sqrt(mean((y_test - gee_rf_predict)^2)))

plot(y_test, rf_predict, xlim = lims, ylim = lims, main = sprintf("R randomForest (RMSE = %.4f)", rf_rmse))
lines(c(-9999, 9999), c(-9999, 9999), col = "red")


sklearn <- reticulate::import("sklearn")
rf_sklearn <- sklearn$ensemble$RandomForestRegressor(n_estimators = as.integer(100), max_features = 0.333333)$fit(x, y)
sk_predict <- rf_sklearn$predict(x_test)
(sk_rmse <- sqrt(mean((y_test - sk_predict)^2)))

plot(y_test, sk_predict, xlim = lims, ylim = lims, main = sprintf("sklearn randomForest (RMSE = %.4f)", sk_rmse))
lines(c(-9999, 9999), c(-9999, 9999), col = "red")

# Convert to gee
geemap <- import("geemap")
trees <- geemap$ml$rf_to_strings(rf_sklearn, feature_names = names(x_test), processes = as.integer(8), output_mode = "REGRESSION")

writeLines(trees, "output.txt")
trees
ee_classifier <- geemap$ml$strings_to_classifier(res)
gee_predict2 <- ee_classifier %>% predict(test_extract)
(gee_rmse2 <- sqrt(mean((y_test - gee_predict2)^2)))

writeLines(trees, "../trees.txt")

library(randomForest)
data(iris)
rf_iris <- randomForest(iris[, 1:3], iris[, 4])
randomForest::getTree(rf_iris, 1)
sum(df$left == 0)
sum(df$left != 0)

res <- build_forest(rf)
strsplit(res[1], "\n")[[1]][1:10]
