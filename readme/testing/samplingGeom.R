# devtools::load_all()
source("readme/testing/load_atl08_data.R")

index <- ANNIndex$new(all_data$longitude, all_data$latitude)

geomS <- geomSampling(0.9999, geom = geom)
stratS <- stratifiedSampling(3, "h_canopy", breaks = c(0, 10, 20, Inf))

sampled <- sample(
    all_data,
        geomSampling(0.9999999, geom = geom) +
        stratifiedSampling(0.99999, "h_canopy", breaks = c(0, 10, 20, Inf)) +
    spacedSampling(3, radius = 0.01, spatialIndex = index) 
)


library(ggplot2)
# install.packages('tidyterra')
library(tidyterra)

ggplot() +
    tidyterra::geom_spatvector(data = geom) +
    # tidyterra::geom_spatvector(data=geom_v, aes(color = as.factor(rowid))) +
    ggplot2::geom_point(data = all_data, aes(longitude, latitude, color = h_canopy)) +
    theme_bw()


ggplot() +
    tidyterra::geom_spatvector(data = geom) +
    # tidyterra::geom_spatvector(data=geom_v, aes(color = as.factor(rowid))) +
    ggplot2::geom_point(data = sampled, aes(longitude, latitude, color = breaks)) +
    theme_bw()
