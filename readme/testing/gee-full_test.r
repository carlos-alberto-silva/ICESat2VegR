# source("C:/Users/caiohamamura/src/r/rICESat2Veg/readme/testing/gee-full_test.r")
library(terra)
library(magrittr)
library(ICESat2VegR)

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
all_data <- NULL
for (granule in granulesName) {
  atl08 <- ICESat2VegR::ATL08_read(file.path("../inst/exdata", granule))
  dt <- ICESat2VegR::ATL08_seg_attributes_dt(atl08, attribute = c("h_canopy"))
  dt <- suppressWarnings(dt[, granule := granule])
  if (is.null(all_data)) {
    all_data <- dt
  } else {
    all_data <- c(dt, all_data)
  }
}


all_data <- all_data %>% ATL08_seg_attributes_dt_clipGeometry(geom)

radius <- 0.01
# Build a new Spatial Index
ann <- ANNIndex$new(all_data$longitude, all_data$latitude)

library(magrittr)

# Get some samples
set.seed(123)
sample1 <- sample(all_data, gridSampling) 
sample2 <- sample(sample1, stratifiedSampling())


library(ggplot2)
library(terra)
library(tidyterra)

ggplot() + 
  geom_spatvector(data = geom, color = "orange", lwd=2, fill = NA) +
  ggplot2::geom_point(data = the_sample, aes(longitude, latitude)) + 
  theme_bw()

