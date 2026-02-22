library(ICESat2VegR)

# Specifying bounding box coordinates
lower_left_lon <- -84.0
lower_left_lat <- 30.0
upper_right_lon <- -82.5
upper_right_lat <- 34.0

# Specifying the date range
daterange <- c("2019-04-12", "2019-05-03")
clip_region <- terra::ext(-83.26, -83.11, 31.95, 32.46)

atl03_granules_local <- ATLAS_dataFinder(
  short_name = "ATL03",
  clip_region$xmin,
  clip_region$ymin,
  clip_region$xmax,
  clip_region$ymax,
  version = "007",
  daterange = daterange,
  persist = FALSE,
  cloud_computing = FALSE
)

atl08_granules_local <- ATLAS_dataFinder(
  short_name = "ATL08",
  lower_left_lon,
  lower_left_lat,
  upper_right_lon,
  upper_right_lat,
  version = "007",
  daterange = daterange,
  persist = TRUE,
  cloud_computing = FALSE
)

# Download granules
ATLAS_dataDownload(atl03_granules_local, outdir)
ATLAS_dataDownload(atl08_granules_local, outdir)


## ATL03
# Read the granules
atl03_files <- list.files(outdir, "ATL03.*h5", full.names = TRUE)
atl03_h5 <- lapply(atl03_files, ATL03_read)

## ATL08
# Read the granules
atl08_files <- list.files(outdir, "ATL08.*h5", full.names = TRUE)
atl08_h5 <- lapply(atl08_files, ATL08_read)


for (ii in c(1,2)) {
    output <- sub(".h5", "_clip.h5", atl03_files[ii])
    atl03_clipped <- ATL03_h5_clipBox(atl03_h5[[ii]], output, clip_obj = clip_region, beam = c("gt1r", "gt2r", "gt3r"), additional_groups=c("METADATA", "ancillary_data", "orbit_info", "quality_assessment"))
    close(atl03_clipped)
}
for (ii in c(1,2)) {
    output <- sub(".h5", "_clip.h5", atl08_files[ii])
    atl08_clipped <- ATL08_h5_clipBox(atl08_h5[[ii]], output, clip_obj = clip_region, beam = c("gt1r", "gt2r", "gt3r"), additional_groups=c("METADATA", "ancillary_data", "orbit_info", "quality_assessment"))
    close(atl08_clipped)
}