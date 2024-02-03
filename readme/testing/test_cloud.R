library(rICESat2Veg)

lower_left_lon <- -85
lower_left_lat <- 30.0
upper_right_lon <- -84.0
upper_right_lat <- 31.0
version <- "006"
short_name <- "ATL08"
persist <- TRUE
daterange <- c("2019-08-19", "2019-08-20")


res <- ATLAS_dataFinder_cloud(
  short_name,
  lower_left_lon,
  lower_left_lat,
  upper_right_lon,
  upper_right_lat,
  version,
  daterange,
  persist
)

f <- ATL03_read(res[1])

test <- ATL03_photons_attributes_dt(f)

if (!auth$authenticated) {
  auth <- ea$login(strategy = "environment")
}

if (!auth$authenticated) {
  auth <- ea$login(strategy = "interactive", persist = persist)
}

granules <- ATL03_read(res[[1]])
ATL03_read("inst")
f <- h5py$File(granules[[1]])
g = f['gt1l']
l = g[['land_segments/latitude']]

atl08_granules <- ATLAS_dataFinder(
  "ATL08", ul_lat, ul_lon, lr_lat, lr_lon, version, daterange
)

output_dir <- tempdir()
setwd(output_dir)

ATLAS_dataDownload(atl03_granules[1], ".")
ATLAS_dataDownload(atl08_granules[1], ".")

atl03_path <- list.files(pattern = "^ATL03")[1]
atl08_path <- list.files(pattern = "^ATL08")[1]

atl03_h5 <- ATL03_read(atl03_path)
atl08_h5 <- ATL08_read(atl08_path)

atl03_atl08_dt <- ATL03_ATL08_photons_attributes_dt_join(atl03_h5, atl08_h5, beam = "gt1l")
plot(atl03_atl08_dt, "h_ph", xlim = c(406000, 410000), ylim = c(180, 220), pch = 20, cex = 0.7)

dt <- data.table::as.data.table(
  na.omit(atl03_atl08_dt)
)[dist_ph_along > 406000 & dist_ph_along < 410000]



atl03_atl08_dt[dist]
library(ggplot2)

ggplot(
  dt,
  aes(
    x = dist_ph_along,
    y = h_ph,
    colour = as.factor(classed_pc_flag)
  )
) +
  geom_point() +
  scale_color_manual(
    values = c(
      "0" = "gray",
      "1" = "brown",
      "2" = "green",
      "3" = "forestgreen"
    )
  )
