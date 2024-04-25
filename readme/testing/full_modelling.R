# source("R/gee-search.R")
# source("R/gee-base.R")
# geemap <- reticulate::import("geemap")

devtools::load_all()
library(terra)
library(magrittr)

geom <- terra::vect("../inst/extdata/all_boundary.shp")
bbox <- terra::ext(geom)

years <- c(2019)
aprilPlaceholder <- "%s-04-01"
mayPlaceholder <- "%s-05-31"

all_granules <- c()
year <- 2019

for (year in years) {
  granules <- ICESat2VegR::ATLAS_dataFinder(
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

strong_beams <- c(
  "gt1r",
  "gt2r",
  "gt3r"
)

target_attributes <- c("h_canopy")

all_dt <- list()
ii <- 1
for (granule in granules) {
  atl08 <- ATL08_read(granule)
  all_dt[[ii]] <- ATL08_seg_attributes_dt(atl08, beam = strong_beams, attribute = target_attributes)[h_canopy < 100]
  ii <- ii + 1
}
dt <- data.table::rbindlist(all_dt)
dt
prepend_class(dt, "icesat2.atl08_dt")

dt2 <- ATL08_seg_attributes_dt_clipGeometry(dt, polygon = geom, split_by = "layer")
nrow(dt2)


#####################
## EARTH ENGINE
#####################
source("readme/testing/02gee_modelling")



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
  as.data.frame(sampled[, .SD, .SDcols = c("beam", "longitude", "latitude", "h_canopy")]),
  geom = c("longitude", "latitude")
)

names(sampled_vect)

ii <- 1

# fullStack = hls
ee_sampled <- extract(fullStack, sampled_vect, 30)

INT_MAX <- 2147483647



#### INPUT ################
x <- ee_sampled
y_name <- "h_canopy"
nTrees <- 100
train_split <- 0.7
###########################

selected_properties <- ee$List(list())
ee_bandNames <- x$first()$propertyNames()
ee_bandNames <- ee_bandNames$remove("system:index")$remove("h_canopy")$remove("beam")
n <- ee_sampled$size()
train_size <- n$multiply(train_split)$round()$int()
validation_size <- n$subtract(train_size)

result <- list(
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
      return(f)
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


pymain <- reticulate::import_main()



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
