# source("R/gee-search.R")
# source("R/gee-base.R")
# geemap <- reticulate::import("geemap")


library(ICESat2VegR)
library(terra)
library(magrittr)

geom <- terra::vect("../output.geojson")
ext <- terra::ext(geom)
center <- as.numeric(terra::geom(terra::centroids(terra::vect(ext)))[, c("x", "y")])
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

bandNames(collection)
pybase <- reticulate::import_builtins()
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



img <- hls
visBand <- "blue"

img <- blue$visualize(
  bands = c(visBand),
  min = min(img[[visBand]]),
  max = max(img[[visBand]])
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





source("R/extract.R")

# Extract
train <- terra::vect("../train2.geojson")
test <- terra::vect("../test.geojson")

train <- train[train$h_canopy <= 50]
test <- test[test$h_canopy <= 50]

train_extract <- extract(fullStack, train, 30)
test_extract <- extract(fullStack, test, 30)
train_df <- ee_to_df(train_extract)
test_df <- ee_to_df(test_extract)
x <- train_df[, -c(1, which(names(train_df) == "h_canopy"))]
y <- train_df[, "h_canopy"]
x_test <- test_df[, -c(1, which(names(train_df) == "h_canopy"))]
y_test <- test_df[, "h_canopy"]

randomForestClassifier <- randomForestRegression(train_extract, property_name = "h_canopy", nTrees = 100, nodesize = 1)
gee_predict <- randomForestClassifier %>% predict(test_extract)
(gee_rmse <- sqrt(mean((y_test - gee_predict)^2)))
lims <- c(0, max(c(y_test, gee_predict)))
plot(y_test, gee_predict, xlim = lims, ylim = lims, main = sprintf("GEE randomForest (RMSE = %.4f)", gee_rmse))
lines(c(-9999, 9999), c(-9999, 9999), col = "red")

library(randomForest)
rf <- randomForest(x, y, ntree = 100, nodesize = 1)
rf_predict <- rf %>% predict(test_df)
rf_rmse <- sqrt(mean((y_test - rf_predict)^2))

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
trees <- geemap$ml$rf_to_strings(rf_sklearn, feature_names = names(x_test), processes = 8, output_mode = "REGRESSION")
ee_classifier <- geemap$ml$strings_to_classifier(trees)
gee_predict2 <- ee_classifier %>% predict(test_extract)
(gee_rmse2 <- sqrt(mean((y_test - gee_predict2)^2)))

writeLines(trees, "../trees.txt")

data(iris)
rf_iris = randomForest(iris[,1:3], iris[,4])
df = data.frame(
  left = rf_iris$forest$leftDaughter[1:87,1],
  right = rf_iris$forest$rightDaughter[1:87,1],
  bestvar = rf_iris$forest$bestvar[1:87,1],
  xbestsplit = rf_iris$forest$ [1:87,1],
  nodestatus = rf_iris$forest$nodestatus[1:87,1],
  nodepred = rf_iris$forest$nodepred[1:87,1]
  )

sum(df$left == 0)
sum(df$left != 0)
