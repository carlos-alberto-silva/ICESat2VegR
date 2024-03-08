# atl08 <- ATL08_read("../inst/extdata/ATL08_20191002161359_00880506_006_02.h5")
# library(ICESat2VegR)
devtools::load_all()
library(data.table)
source("readme/testing/load_atl08_data.R")

selected_name <- all_data[, .N, by=granule][which.max(N), granule]
selected_name <- gsub(".h5", "", selected_name)

bounds <- terra::vect("Z:/01_Projects/04_NASA_ICESat2/00_data/01_SHPs/all_boundary.shp")
ext <- terra::ext(bounds)

granules <- ICESat2_finder(
  short_name = "ATL03",
  lower_left_lon = ext$xmin,
  upper_right_lon = ext$xmax,
  lower_left_lat = ext$ymin,
  upper_right_lat = ext$ymax,
  daterange = c(fall_post_michael_date_start, fall_post_michael_date_end)
)

ATL03_name <- gsub("ATL08", "ATL03", selected_name)
selected_granule <- grep(ATL03_name, granules, value = T)
ICESat2_download(selected_granule, "../inst/extdata/")

granule_id <- gsub("ATL03_", "", ATL03_name)
atl03 <- ATL03_read(gettextf("../inst/extdata/ATL03_%s.h5", granule_id))
atl08 <- ATL08_read(gettextf("../inst/extdata/ATL08_%s.h5", granule_id))


output <- tempfile(fileext = ".h5")
out <- clip(atl03, output, bounds)
atl03_dt <- ATL03_photons_attributes_dt(out)
atl03_dt

output2 <- tempfile(fileext = ".h5")
out2 <- clip(atl08, output2, bounds)
atl08_dt <- ATL08_seg_attributes_dt(out2)

library(tidyterra)
library(ggplot2)

ggplot() + 
  geom_spatvector(data = bounds) +
  geom_point(data = atl03_dt, aes(lon_ph, lat_ph))


ggplot() + 
  geom_spatvector(data = bounds) +
  geom_point(data = atl08_dt, aes(longitude, latitude))

sum(out[["gt1l/geolocation/segment_ph_cnt"]][])

atl03_08_dt <- ATL03_ATL08_photons_attributes_dt_join(out, out2)
atl03_08_dt <- na.omit(atl03_08_dt)
dt <- atl03_08_dt[ph_h < 100 & beam == "gt1r"]
dt <- dt[, classed_pc_flag := factor(classed_pc_flag)]
library(ggplot2)
bbox <- ext
# output = tempfile(fileext=".h5")




{
  min_x <- 0
  max_x <- min_x + 10000
  windows()
  ggplot() +
    geom_point(data = dt[dist_ph_along > min_x & dist_ph_along <= max_x], aes(dist_ph_along, ph_h, color = classed_pc_flag)) +
    ylim(c(0, 25)) +
    xlim(min_x, max_x) +
    scale_color_manual(
      values = c("gray", "orange", "forestgreen", "green"),
      labels = c("Noise", "Ground", "Vegetation", "High Vegetation")
    ) +
    labs(
      color = "Photon classification",
      title = gettextf("Granule: %s", ATL03_name)
    ) +
    theme_bw()
}


