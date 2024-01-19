# install.packages("lidR")

require(lidR)
library(rICESat2Veg)
atl08_path <- "inst\\exdata\\ATL08_20220401221822_01501506_005_01.h5"
atl08_h5 <- ATL08_read(atl08_path)
atl03_path <- "inst\\exdata\\ATL03_20220401221822_01501506_005_01.h5"
atl03_h5 <- ATL03_read(atl03_path = atl03_path)

## join
atl08_dt <- ATL08_seg_attributes_dt(atl08_h5, beam = "gt1l")



output <- "../atl08output.laz"
source('R/utmTools.r')

ATL08_seg_attributes_dt_LAS(atl08_dt, output)

# Record an EPSG code
# epsg(header) <- 32618
