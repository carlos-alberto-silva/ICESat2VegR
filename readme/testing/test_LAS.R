# install.packages("lidR")

require(lidR)
library(ICESat2VegR)


atl08_path <- "../inst\\exdata\\ATL08_20220826063041_09981605_005_01.h5"
# atl03_path <- "inst\\exdata\\ATL03_20220401221822_01501506_005_01.h5"
atl08_h5 <- ATL08_read(atl08_path)
# atl03_h5 <- ATL03_read(atl03_path = atl03_path)

ds <- atl08_h5[["gt1l/land_segments/latitude"]]
ds[]


atl03_atl08_dt <- ATL03_ATL08_photons_attributes_dt_join(atl03_h5, atl08_h5, beam = c("gt1l"))

segment_length <- 30
output <- NA
overwrite <- FALSE


centroid <- "mean"
atl03_atl08_seg_dt <- ATL03_ATL08_segment_create(atl03_atl08_dt, centroid = centroid, segment_length = segment_length, output = output)

output <- "../output.gpkg"
centroid <- "midpoint"
atl03_atl08_seg_dt <- ATL03_ATL08_segment_create(atl03_atl08_dt, centroid = centroid, segment_length = segment_length, output = output)


smoothing_window <- NA
smoothing_func <- median
library(signal)
interpolation_func <- signal::pchip



diff <- max_x - min_x

x <- seq(min_x, max_x, 0.01)


l <- ATL03_ATL08_photons_fitground_seg_dt(
  atl03_atl08_seg_dt,
  smoothing_window,
  smoothing_func,
  interpolation_func,
  xi = atl03_atl08_seg_dt$dist_ph_along
)

l <- ATL03_ATL08_photons_dt_height_normalize(
  atl03_atl08_seg_dt,
  smoothing_window,
  smoothing_func,
  interpolation_func,
  xout_parameter_name = "xi"
)

min_x <- 1850
max_x <- 1900

with(
  atl03_atl08_seg_dt[
    classed_pc_flag >= 1 &
      dist_ph_along >= min_x &
      dist_ph_along < max_x
  ],
  plot(
    dist_ph_along, ph_h, 
    xlim = c(min_x, max_x), 
    col = c("orange", "forestgreen", "green")[classed_pc_flag])
)
which.min(l[
    classed_pc_flag >= 1
  , ph_h])

mid_x = l[classed_pc_flag >= 1][7837932, dist_ph_along]
min_x = mid_x - 100
max_x = mid_x + 100

windows()
with(
  l[
      dist_ph_along >= min_x &
      dist_ph_along < max_x
  ],
  plot(
    dist_ph_along, ph_h, 
    xlim = c(min_x, max_x), 
    col = c("gray","orange", "forestgreen", "green")[classed_pc_flag+1])
)

windows()
with(
  atl03_atl08_dt[
      dist_ph_along >= min_x &
      dist_ph_along < max_x
  ],
  plot(
    dist_ph_along, ph_h, 
    xlim = c(min_x, max_x), 
    col = c("gray","orange", "forestgreen", "green")[classed_pc_flag+1])
)

windows()
with(
  atl03_atl08_dt[
      dist_ph_along >= min_x &
      dist_ph_along < max_x
  ],
  plot(
    dist_ph_along, h_ph, 
    xlim = c(min_x, max_x), 
    col = c("gray","orange", "forestgreen", "green")[classed_pc_flag+1])
)


plot(x, l)
atl03_atl08_seg_dt[
  dist_ph_along >= min_x &
    dist_ph_along < max_x &
    classed_pc_flag == 1,
  points(dist_ph_along, h_ph)
]

# ## join
# atl03_atl08_dt <- ATL03_ATL08_photons_attributes_dt_join(atl03_h5, atl08_h5, beam = "gt1l")
# # atl08_dt <- ATL08_seg_attributes_dt(atl08_h5, beam = "gt1l")

# atl08_dt <- ATL08_seg_attributes_dt(atl08_h5)


# output <- "../atl08output.laz"
# source('R/utmTools.r')

# ATL08_seg_attributes_dt_LAS(atl08_dt, output)

# # Record an EPSG code
# # epsg(header) <- 32618

# ground_track
# atl03_dt <- ATL03_photons_attributes_dt(atl03_h5, beam = "gt1l")

# output <- "../rgt.gpkg"
# gt <- ATL03_rgt_extract(atl03_h5, output = output)

# las <- ATL03_ATL08_photons_attributes_dt_LAS(atl03_atl08_seg_dt, output = "../output.las")
