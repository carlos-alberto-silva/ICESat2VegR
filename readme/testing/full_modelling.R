# source("R/gee-search.R")
# source("R/gee-base.R")
# geemap <- reticulate::import("geemap")

devtools::load_all()
library(terra)
library(magrittr)

geom <- terra::vect("Z:/01_Projects/04_NASA_ICESat2/00_data/01_SHPs/all_boundary.shp")
bbox <- terra::ext(geom)

years <- c(2019:2022)
aprilPlaceholder <- "%s-04-01"
mayPlaceholder <- "%s-05-31"

all_granules <- c()

for (year in years) {
  granules <- ICESat2VegR::ICESat2_finder(
    short_name = "ATL08",
    version = "006",
    daterange = c(gettextf(aprilPlaceholder, year), gettextf(mayPlaceholder, year)),
    lower_left_lon = bbox$xmin,
    lower_left_lat = bbox$ymin,
    upper_right_lon = bbox$xmax,
    upper_right_lat = bbox$ymax
  )
  all_granules <- c(all_granules, granules)
}


ICESat2_download(all_granules, "../inst/extdata")


current_year <- 2019
granules <- list.files("../inst/extdata", "2019", full.names = TRUE)

power_beams <- c(
  "gt1r",
  "gt2r",
  "gt3r"
)

target_attributes <- c("h_canopy")

all_dt <- list()
ii <- 1
for (granule in granules) {
  atl08 <- ATL08_read(granule)
  all_dt[[ii]] <- ATL08_seg_attributes_dt(atl08, beam = power_beams, attribute = target_attributes)[h_canopy < 100]
  ii <- ii + 1
}
dt <- data.table::rbindlist(all_dt)
dt
prepend_class(dt, "icesat2.atl08_dt")

dt2 <- ATL08_seg_attributes_dt_clipGeometry(dt, polygon = geom, split_by = "layer")
nrow(dt2)


aoi <- ee$Geometry$BBox(
  west = bbox$xmin,
  south = bbox$ymin,
  east = bbox$xmax,
  north = bbox$ymax
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
  median() $
  clip(aoi)



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
result <- search_datasets("nasa", "dem")
catalog_id <- get_catalog_id(result)

elevation <- ee$Image(catalog_id)

the_slope <- as.integer(slope(as.integer(elevation)) * 1000)
the_aspect <- elevation %>% aspect()

# stackDem
stackDem <- c(elevation, the_slope, the_aspect)$clip(aoi)


#############################
## SENTINEL 1C
#############################

s1c = ee$ImageCollection("COPERNICUS/S1_GRD")$
  filterBounds(aoi)$
  filterDate("2020-02-01", "2020-05-31")$
  filter(ee$Filter$listContains('transmitterReceiverPolarisation', 'VV'))$
  filter(ee$Filter$listContains('transmitterReceiverPolarisation', 'VH'))$
  filter(ee$Filter$eq('instrumentMode', 'IW'))

s1c = s1c$sort('system:time_start', FALSE)
s1c = s1c$reduce(ee$Reducer$firstNonNull())
s1c = s1c[['VV_first', 'VH_first']]
names(s1c) = c('vv', 'vh')
s1c = as.integer(s1c * 100)

# Radar vegetation index
vv = s1c[["vv"]]
vh = s1c[["vh"]]

# RVI
s1c[["rvi"]] = as.integer(sqrt(vv/(vv+vh))*(vv/vh) * 10000)

# COPOLs
s1c[["copol"]] = as.integer((vv / vh) * 10000)
s1c[["copol2"]] = as.integer(((vv-vh)/(vv+vh)) * 10000)
s1c[["copol3"]] = as.integer(vh / vv) * 10000)



##################
## FINAL STACK
##################
fullStack <- c(hls, s1c, stackDem)

kernel <- ee$Kernel$fixed(3, 3, list(c(1, 1, 1), c(1, 1, 1), c(1, 1, 1)), 3, 3, FALSE)
reducerNames <- c("mean", "min", "max", "stdDev")
for (reducerName in reducerNames) {
  reducer <- ee$Reducer[[reducerName]]()
  fullStack[[""]] <- hls$reduceNeighborhood(reducer, kernel)
  fullStack[[""]] <- stackDem$reduceNeighborhood(reducer, kernel)
}
fullStack


# source("R/gee-base.R")

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
center <- terra::geom(terra::centroids(geom))


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
  leaflet::setView(lng = center[, "x"][[1]], lat = center[, "y"][[1]], zoom = 11)


leaflet_map

# Extract
# all_data_vect <- terra::vect(
#   as.data.frame(dt2),
#   geom = c("longitude", "latitude")
# )

# terra::writeVector(all_data_vect, "../inst/extdata/all_data.geojson", filetype = "geojson")


# all_extract <- list()

# for (ii in seq(1, nrow(all_data_vect), 1000)) {
#   message(gettextf("Extracting rows %d-%d (total %d)", ii, ii + 999, nrow(all_data_vect)))
#   to_extract <- extract(fullStack, all_data_vect[ii:(ii + 999)], 30)
#   all_extract[[""]] <- ee_to_dt(to_extract)
# }

# extract_dt[, "system:index" := gsub("_0", "", `system:index`)]

# head(extract_dt[["system:index"]])
# numeric_dt <- extract_dt[, lapply(.SD, as.numeric), .SDcols = -"beam"]
# numeric_dt[, beam := extract_dt$beam]

# extract_dt <- data.table::rbindlist(all_extract)
# saveRDS(numeric_dt, "../inst/extdata/all_extracted.rds")
# write.csv(numeric_dt, "../inst/extdata/all_extracted.csv")


## LOAD SAVED DATA
set.seed(47289143)
degree_to_meter_factor <- 111139
sampled <- sample(dt2, method = spacedSampling(200, radius = 30 / degree_to_meter_factor))
nrow(sampled)


sampled_vect <- terra::vect(
  as.data.frame(sampled[,.SD,.SDcols=c("beam","longitude","latitude","h_canopy")]),
  geom = c("longitude", "latitude")
)

names(sampled_vect)

ii = 1

# fullStack = hls
ee_sampled <- extract(fullStack, sampled_vect, 30)

INT_MAX = 2147483647



#### INPUT ################ 
x = ee_sampled
y_name = "h_canopy"
nTrees = 100
train_split = 0.7
###########################

selected_properties = ee$List(list())
ee_bandNames <- x$first()$propertyNames()
ee_bandNames <- ee_bandNames$remove("system:index")$remove("h_canopy")$remove("beam")
n <- ee_sampled$size()
train_size <- n$multiply(train_split)$round()$int()
validation_size = n$subtract(train_size)

result = list(
    property = ee$List(list()),
    rmse = ee$List(list())
  )

set.seed(47289143)
for (band in ee_bandNames$getInfo()) {
  current <- selected_properties$add(band)
  
  rmseList <- ee$List(list())
  message(gettextf("Testing %s", current$getInfo()), appendLF = TRUE)
  for (i in 1:100) {
    message(gettextf("\rBoot %d/%d", i, 100), appendLF = FALSE)
    x <- x$randomColumn(seed = floor(runif(1) * INT_MAX))
    train_sample <- x$limit(train_size, "random")$select(current$add(y_name))
    validation_sample <- x$limit(validation_size, "random", ascending = FALSE)$select(current$add(y_name))
    
    
    randomForestClassifier <- randomForestRegression(train_sample, property_name = y_name, train_properties = current, nTrees = nTrees)
    classification <- validation_sample$classify(randomForestClassifier)
    classification2 <- classification$map(function(f) {
      sqerror <- f$getNumber("h_canopy")$subtract(f$getNumber("classification"))$pow(2)
      f <- f$set(list(sqerror = sqerror))
      return (f)
    })
    rmse <- classification2$aggregate_mean("sqerror")$sqrt()
    rmseList <- rmseList$add(rmse)
  }
  message(appendLF = TRUE)
  
  meanRmse <- rmseList$reduce(ee$Reducer$mean())
  result$rmse <- result$rmse$add(meanRmse)
  result$property <- result$property$add(band)
}
minRmse <- result$rmse$reduce(ee$Reducer$min())
minIdx <- result$rmse$indexOf(minRmse)
ee$Dictionary(result)$getInfo()
result$property$getString(minIdx)$getInfo()


pymain = reticulate::import_main()



test_sample <- ee_sampled$limit(validation_size, "random", ascending = FALSE)$select(current)
gee_predict <- randomForestClassifier %>% predict(test_sample)



StatModel(y_test, gee_predict2, xlim = c(0, 30), ylim = c(0, 30))

ee_sampled


# library(data.table)
# numeric_dt <- data.table::as.data.table(readRDS("../inst/extdata/all_extracted.rds"))
# head(numeric_dt)
# range(numeric_dt)
# x_columns <- names(numeric_dt[, .SD, .SDcols = -c("nid", "ids", "system:index", "h_canopy")])
# x <- numeric_dt[, .SD, .SDcols = x_columns]
# y <- numeric_dt$h_canopy

# # Remove columns with no data (all 0)
# all0 <- x[, lapply(.SD, function(x) all(x == 0))]
# remove_cols <- names(all0)[as.logical(all0)]
# x[, eval(remove_cols) := NULL]



# (gee_rmse <- sqrt(mean((y_test - gee_predict)^2)))
# lims <- c(0, max(c(y_test, gee_predict)))
# plot(y_test, gee_predict, xlim = lims, ylim = lims, main = sprintf("GEE randomForest (RMSE = %.4f)", gee_rmse))
# lines(c(-9999, 9999), c(-9999, 9999), col = "red")



library(randomForest)


rf <- randomForest::randomForest(x, y, ntree = 50)
rf
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
