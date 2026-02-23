# ── Machine Learning: Predict and Rasterize ───────────────────────────────────
library(ICESat2VegR)
library(randomForest)
library(purrr)

setwd("~/ic2_workshop_data")

# Sample training data
h_canopy_train <- c(
  13.9, 3.1, 2.2, 4.6, 21.6, 7.2, 5, 7.7, 0.8, 9.7,
  11, 11.3, 15.5, 5.1, 10.4, 0.6, 14.6, 13.3, 9.8, 14.7
)
agbd_train <- c(
  144.8, 27.5, 51.6, 60.5, 232.3, 102.8, 33.1, 91.3, 23,
  120.1, 125.7, 127.2, 147.4, 48.8, 103.3, 55.9, 181.8, 139.9, 120.1, 162.8
)

set.seed(42)
model <- randomForest(data.frame(h_canopy = h_canopy_train), agbd_train)

# Predict for all ATL08 data and write to HDF5
atl08_files <- list.files(pattern = "ATL08.*clip.h5$")
atl08_h5 <- atl08_files |> map(ATL08_read)


out_h5 <- tempfile(fileext = ".h5")
predicted_h5 <- NULL
# Now we will predict the data for the entire ATL08 dataset,
# this can be done in batches if the dataset is too large to fit in memory
# and the results will be appended to the same HDF5 file.
for (x in atl08_h5) {
  seg_dt <- ATL08_seg_attributes_dt(x, attributes = "h_canopy")
  seg_dt[h_canopy > 100, h_canopy := NA_real_]
  seg_dt <- na.omit(seg_dt)
  predicted_h5 <- predict_h5(model, seg_dt, out_h5)
}

# Rasterize predictions
output_raster <- tempfile(fileext = ".tif")
x <- predicted_h5[["longitude"]][]
y <- predicted_h5[["latitude"]][]
bbox <- terra::ext(min(x), max(x), min(y), max(y))


predicted_h5_dt <- data.table::data.table(
  longitude = predicted_h5[["longitude"]][],
  latitude = predicted_h5[["latitude"]][],
  predicted_agbd = predicted_h5[["prediction"]][]
)

terra::vect(predicted_h5_dt, geom = c("longitude", "latitude"), crs = 4326) |>
  terra::writeVector("predicted_agbd.gpkg")

# Rasterize the predicted HDF5 data to a GeoTIFF file
# The points that fall into the same pixel will create
# a tif file with 6 bands: count, sum of squares,
# min, max and standard deviation.
# The sum of squares is used to calculate the standard deviation
# using an online algorithm that does not require storing all the
# values in memory.
rasterize_h5(predicted_h5, "output.tif", bbox = bbox, res = 0.005)

# Plot
# Open the raster by file path
forest_height_palette <- c("#ffffff", "#4d994d", "#004d00")

# Open band 2 only (mean AGBD)
stars_rast <- stars::read_stars(output_raster, RasterIO = list(bands = 2))
res_map <- mapview::mapview(
  stars_rast,
  layer.name = "AGBD mean",
  col.regions = forest_height_palette,
  na.alpha = 0.1,
  map = leaflet::leaflet() |> leaflet::addProviderTiles("Esri.WorldImagery")
)

res_map@map
