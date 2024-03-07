devtools::load_all()
source("readme/testing/load_atl08_data.R")

sample_rast <- terra::rast("Z:/01_Projects/04_NASA_ICESat2/00_data/02_TIFs/split_sample.tif")
sampled <- sample(all_data, rasterSampling(10, sample_rast))

library(ggplot2)
library(tidyterra)

ggplot() +
    tidyterra::geom_spatraster(data = sample_rast) +
    # tidyterra::geom_spatvector(data=geom_v, aes(color = as.factor(rowid))) +
    ggplot2::geom_point(data = sampled, aes(longitude, latitude, color = group)) +
    theme_bw()
