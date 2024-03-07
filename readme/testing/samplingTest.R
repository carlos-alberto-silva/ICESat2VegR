# source("C:/Users/caiohamamura/src/r/rICESat2Veg/readme/testing/gee-full_test.r")
# devtools::load_all()
source("readme/testing/load_atl08_data.R")

radius <- 0.01

# Build a new Spatial Index
ann <- ANNIndex$new(all_data$longitude, all_data$latitude)

# devtools::load_all()
gridSamp <- gridSampling(size = 0.999999, grid_size = 1)
stratSamp <- stratifiedSampling(size = 2, variable = "h_canopy", breaks = c(0, 10, 20, 9999))
spacedSamp
combinedSamp <- gridSamp + stratSamp
# chainSampling = stratSamp
?graphics::hist
(sampled <- sample(all_data, method = combinedSamp))
unique(sampled$breaks)


library(ggplot2)
library(terra)
library(tidyterra)

ggplot() +
  geom_spatvector(data = geom, color = "orange", fill = NA) +
  geom_point(data = sampled, aes(longitude, latitude, color = h_canopy), cex = 3) +
  scale_x_continuous(breaks = seq(min(all_data$longitude), max(all_data$longitude), 1)) +
  scale_y_continuous(breaks = seq(min(all_data$latitude), max(all_data$latitude), 1)) +
  scale_color_gradient(low = "gray", high = "forestgreen") +
  ggplot2::theme_bw()
