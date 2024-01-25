# install.packages("lidR")

require(lidR)
library(rICESat2Veg)


atl08_path <- "inst\\exdata\\ATL08_20220401221822_01501506_005_01.h5"
atl03_path <- "inst\\exdata\\ATL03_20220401221822_01501506_005_01.h5"
atl08_h5 <- ATL08_read(atl08_path)
atl03_h5 <- ATL03_read(atl03_path = atl03_path)

output <- "../output.gpkg"


dt <- ATL03_ATL08_photons_attributes_dt_join(atl03_h5, atl08_h5, beam = c("gt1r"))

segment_length <- 30

ATL03_ATL08_segment_create(dt, segment_length = 30, output = "../segments.gpkg")

# ## join
# atl03_atl08_dt <- ATL03_ATL08_photons_attributes_dt_join(atl03_h5, atl08_h5, beam = "gt1l")
# # atl08_dt <- ATL08_seg_attributes_dt(atl08_h5, beam = "gt1l")

# atl08_dt <- ATL08_seg_attributes_dt(atl08_h5)


# output <- "../atl08output.laz"
# source('R/utmTools.r')

# ATL08_seg_attributes_dt_LAS(atl08_dt, output)

# # Record an EPSG code
# # epsg(header) <- 32618
