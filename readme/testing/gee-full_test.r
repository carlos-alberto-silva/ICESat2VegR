# source("C:/Users/caiohamamura/src/r/rICESat2Veg/readme/testing/gee-full_test.r") 
library(ICESat2VegR)
library(terra)
library(magrittr)

geom <- terra::vect("Z:/01_Projects/04_NASA_ICESat2/00_data/01_SHPs/all_boundary.shp")
ext <- terra::ext(geom)



hurricane_michael_date <- "2018-10-10"
fall_pre_michael_date_start <- "2018-09-23"
fall_pre_michael_date_end <- "2018-10-09"

fall_post_michael_date_start <- "2019-10-01"
fall_post_michael_date_end <- "2019-10-31"

granules <- ICESat2VegR::ICESat2_finder(
  short_name = "ATL08",
  version = "006",
  daterange = c(fall_post_michael_date_start, fall_post_michael_date_end),
  lower_left_lon = ext$xmin,
  lower_left_lat = ext$ymin,
  upper_right_lon = ext$xmax,
  upper_right_lat = ext$ymax
)
granules

ICESat2VegR::ICESat2_download(granules, outdir = "../inst/exdata/")

granulesName <- basename(granules)
all_data <- list()
for (granule in granulesName) {
  atl08 <- ICESat2VegR::ATL08_read(file.path("../inst/exdata", granule))
  dt <- ATL08_seg_attributes_dt(atl08, attribute = c("h_canopy"))
  dt <- suppressWarnings(dt[, granule := granule])
  all_data[[granule]] <- dt
}
dt <- data.table::rbindlist(all_data)

ICESat2VegR::prepend_class(dt, "icesat2.atl08_dt")

dt2 <- ATL08_seg_attributes_dt_clipGeometry(dt, geom)

library(Rcpp)
setwd("C:/Users/caiohamamura/src/r/rICESat2Veg/")
gridCpp <- Rcpp::sourceCpp("src/GridIndex.cpp")
Rcpp::loadModule("grid_index_module")
Rcpp::sourceCpp("src/findRadius.cpp")
  
grid <- GridIndex$new(dt2$longitude, dt2$latitude, 0.01)
dt2[,"I" := .I]
(idx = sample(dt2$I, 1))
pt = dt2[idx]

microbenchmark::microbenchmark(grid$searchFixedRadius(pt$longitude, pt$latitude, 0.1), times=1000)

microbenchmark::microbenchmark(findRadius3(dt2$longitude, dt2$latitude, 0.1, 100), times = 100)
findRadius(dt2$longitude, dt2$latitude, 0.1, 100, grid$getPointer())



library(Rcpp)


queryX = 2.5
queryY = 6.5
radius = 1.0

grid$searchFixedRadius(pt$longitude, pt$latitude, 0.01)
grid$searchFixedRadius(queryX, queryY, 2.5)

