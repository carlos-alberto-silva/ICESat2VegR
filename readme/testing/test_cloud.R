library(ICESat2VegR)

lower_left_lon <- -85
lower_left_lat <- 30.0
upper_right_lon <- -84.0
upper_right_lat <- 31.0
version <- "006"
short_name <- "ATL08"
persist <- TRUE
daterange <- c("2019-08-19", "2019-08-20")

# devtools::load_all()
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

atl08_h5 <- ATL08_read(res[1])


# # ATL08_photon.var.map
# ATL08_photon.var.map <- list()
# ATL08_photon.var.map[["ph_segment_id"]] <- "ph_segment_id"
# ATL08_photon.var.map[["classed_pc_indx"]] <- "classed_pc_indx"
# ATL08_photon.var.map[["classed_pc_flag"]] <- "classed_pc_flag"
# ATL08_photon.var.map[["ph_h"]] <- "ph_h"
# ATL08_photon.var.map[["d_flag"]] <- "d_flag"
# ATL08_photon.var.map[["delta_time"]] <- "delta_time"

# beam = c("gt1l", "gt1r", "gt2l", "gt2r", "gt3l", "gt3r")
# photon_attribute = c("ph_segment_id", "classed_pc_indx", "classed_pc_flag", "ph_h", "d_flag", "delta_time")
# # Check beams to select
#   groups_id <- getBeams(atl08_h5)


#   beam <- intersect(groups_id, beam)
#   photon.dt <- list()

#   pb <- utils::txtProgressBar(min = 0, max = length(beam), style = 3)

#   i_s <- 0

#     i = beam[1]
#       i_s <- i_s + 1

#       atl08_h5v2_i <- atl08_h5[[paste0(i, "/signal_photons")]]

#       m <- data.table::data.table()

#         col = photon_attribute[1]
#         metric_address <- ATL08_photon.var.map[[col]]

#         is.null(metric_address)
#           h5exists(atl08_h5v2_i, col)
#             metric_address <- col
          
#         base_addr <- gsub("^(.*)/.*", "\\1", metric_address)
#         h5exists(atl08_h5v2_i, base_addr) && h5exists(atl08_h5v2_i, metric_address)
#           m[, eval(col) := atl08_h5v2_i[[metric_address]][]]


#       photon_dt <- data.table::rbindlist(list(photon.dt, m), fill = TRUE)
#       utils::setTxtProgressBar(pb, i_s)

#   setattr(photon_dt, "class", c("icesat2.atl08_dt", "data.table", "data.frame"))

#   close(pb)






test <- ATL08_photons_attributes_dt(atl08_h5)





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
