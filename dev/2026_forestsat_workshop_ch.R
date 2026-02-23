# ICESat2VegR Workflow Script
# Authors: Carlos Alberto Silva and Caio Hamamura
# Processes NASA ICESat-2 ATL03 and ATL08 data for land and vegetation applications

# ── Installation ──────────────────────────────────────────────────────────────
pak::pkg_install("carlos-alberto-silva/ICESat2VegR", dependencies = TRUE)

# -- Load
library(ICESat2VegR)

# ── Authentication ────────────────────────────────────────────────────────────
# reticulate::py_require("earthengine-api")
# reticulate::py_module_available("earthaccess")
# Option 1: Set env vars EARTHDATA_USERNAME and EARTHDATA_PASSWORD before opening R
# Option 2:

earthdata_login(output_dir = "~")

# Google Earth Engine (if needed)
# Create a project with Earth Engine
ee_initialize("forestsat-workshop-2026")


# ── Download Example Dataset ──────────────────────────────────────────────────

outdir <- tempdir()
dir.create(outdir, showWarnings = FALSE)

url <- "https://github.com/carlos-alberto-silva/ICESat2VegR/releases/download/example_datasets/Study_Site.zip"
zip_file <- file.path(outdir, "Study_Site.zip")
if (!file.exists(zip_file)) download.file(url, zip_file, mode = "wb")

if (!file.exists(file.path(outdir, "ATL03_20190413122158_02330302_006_02_clip.h5"))) {
  unzip(zip_file, exdir = outdir)
}


# ── Search Parameters ─────────────────────────────────────────────────────────
setwd(outdir)
vector_data <- terra::vect('example_aoi.gpkg')
bounding_box <- terra::ext(vector_data)
bounding_box

lower_left_lon  <- bounding_box$xmin
lower_left_lat  <- bounding_box$ymin
upper_right_lon <- bounding_box$xmax
upper_right_lat <-bounding_box$ymax
daterange <- c("2019-04-12", "2019-05-03")


# ── Find & Download Granules (Local Workflow) ─────────────────────────────────

atl03_granules_local <- ATLAS_dataFinder(
  short_name      = "ATL03",
  lower_left_lon, lower_left_lat, upper_right_lon, upper_right_lat,
  version         = "007",
  daterange       = daterange,
  persist         = FALSE,
  cloud_computing = FALSE
)
atl03_granules_local

atl08_granules_local <- ATLAS_dataFinder(
  short_name      = "ATL08",
  lower_left_lon, lower_left_lat, upper_right_lon, upper_right_lat,
  version         = "007",
  daterange       = daterange,
  persist         = TRUE,
  cloud_computing = FALSE
)
atl08_granules_local


# ATLAS_dataDownload(atl03_granules_local, outdir)
# ATLAS_dataDownload(atl08_granules_local, outdir)


# ── Open HDF5 Files ───────────────────────────────────────────────────────────

atl03_files <- list.files(outdir, "ATL03.*h5", full.names = TRUE)
atl03_h5    <- lapply(atl03_files, ATL03_read)

atl08_files <- list.files(outdir, "ATL08.*h5", full.names = TRUE)
atl08_h5    <- lapply(atl08_files, ATL08_read)

# Alternative pipe operator
atl03_h5 |> lapply(close)
atl08_h5 |> lapply(close)

# Alternative purrr
library(purrr)
atl03_h5 <- atl03_files |> map(ATL03_read)
atl08_h5 <- atl08_files |> map(ATL08_read)



# ── Extract ATL03 Photon Attributes ───────────────────────────────────────────
atl03_photons_dt <- atl03_h5 |> map(ATL03_photons_attributes_dt)

atl03_photons_dt[[1]][,name := basename(atl03_files[[1]])]
atl03_photons_dt[[2]][,name := basename(atl03_files[[2]])]


# Quick plot of a subset
atl03_photons_dt[[2]][
  beam == "gt1r" & dist_ph_along < 10000,
  plot(dist_ph_along, h_ph, xlab = "dist_ph_along", ylab = "Elevation (m)", pch = 16, cex = 0.2)
]

atl03_photons_dt <- rbindlist2(atl03_photons_dt)

atl03_photons_dt[
  beam == "gt1r" & dist_ph_along < 10000 & name == "ATL03_20190413122158_02330302_007_01_clip.h5",
  plot(dist_ph_along, h_ph, xlab = "dist_ph_along", ylab = "Elevation (m)", pch = 16, cex = 0.2)
]

# ── Extract ATL03 Segment Metadata ────────────────────────────────────────────

# In a single move
atl03_seg_dt <- atl03_h5 |>
  map(
    ATL03_seg_metadata_dt,
    attributes = c("delta_time", "solar_elevation", "pitch", "h_ph", "ref_elev")
    ) |> rbindlist2()

atl03_seg_dt <- atl03_seg_dt[h_ph < 20000]  # remove outliers


# ── Extract ATL08 Segment Attributes ─────────────────────────────────────────

atl08_seg_dt <- atl08_h5 |>
  map(
    ATL08_seg_attributes_dt,
    attributes = c("h_canopy", "h_te_mean", "terrain_slope", "canopy_openness", "night_flag")
  ) |> rbindlist2()

atl08_seg_dt <- atl08_seg_dt[h_canopy < 100 & h_te_mean < 20000]

# Histograms
layout(t(1:2))
hist(atl03_seg_dt$h_ph,       col = "#bd8421", xlab = "Elevation (m)", main = "ATL03 h_ph")
hist(atl08_seg_dt$h_canopy,   col = "green",   xlab = "Height (m)",    main = "ATL08 h_canopy")
layout(1)

# ── Export to Vector ──────────────────────────────────────────────────────────

atl03_seg_vect <- to_vect(atl03_seg_dt)
atl08_seg_vect <- to_vect(atl08_seg_dt)

# Palette function
greenYellowRed <- function(n) {
  grDevices::hcl.colors(n, "RdYlGn")
}

# Plot with mapview

map_vect <- mapview::mapView(
  atl08_seg_vect,
  layer.name = "h_canopy",
  zcol = "h_canopy",
  col.regions = greenYellowRed,
  map.types = c("Esri.WorldImagery")
)

map_vect

# Write the outputs to gpkg
terra::writeVector(atl03_seg_vect, file.path(outdir, "atl03_seg.gpkg"))
terra::writeVector(atl08_seg_vect, file.path(outdir, "atl08_seg.gpkg"))


# ── Rasterize ATL08 Attributes ────────────────────────────────────────────────

# Single statistic
max_h_canopy <- ATL08_seg_attributes_dt_gridStat(atl08_seg_dt,
  func = max(h_canopy), res = 0.01)


plot(max_h_canopy, main="max h_canopy")

# Multiple statistics
multiple_attributes <- ATL08_seg_attributes_dt_gridStat(atl08_seg_dt, func = list(
  max_h_canopy         = max(h_canopy),
  min_h_canopy         = min(h_canopy),
  mean_canopy_openness = mean(canopy_openness),
  mean_h_te_mean       = mean(h_te_mean)
), res = 0.01)

plot(multiple_attributes)

redYellowGreen <- function(n) grDevices::hcl.colors(n, "RdYlGn")

# With mapview
map_vect_openness <- mapview::mapView(
  atl08_seg_vect,
  zcol = "canopy_openness",
  layer.name = "canopy_openness",
  col.regions = redYellowGreen,
  map.types = c("Esri.WorldImagery")
)

blueYellowRed <- function(n) grDevices::hcl.colors(n, "RdYlBu", rev = TRUE)

map_vect_terrain <- mapview::mapView(
  atl08_seg_vect,
  zcol = "h_te_mean",
  layer.name = "h_te_mean",
  col.regions = blueYellowRed,
  map.types = c("Esri.WorldImagery")
)


m1 <- mapview::mapView(multiple_attributes[[1]], layer.name = "Max h_canopy", map = map_vect, col.regions = redYellowGreen)
m2 <- mapview::mapView(multiple_attributes[[2]], layer.name = "Min h_canopy", map = map_vect, col.regions = redYellowGreen)
m3 <- mapview::mapView(multiple_attributes[[3]], layer.name = "Mean canopy openness", map = map_vect_openness, col.regions = redYellowGreen)
m4 <- mapview::mapView(multiple_attributes[[4]], layer.name = "Mean h_te_mean", col.regions = blueYellowRed, map = map_vect_terrain)

leafsync::sync(m1, m2, m3, m4)

# ── Clipping ──────────────────────────────────────────────────────────────────
# Define bbox
aoi <- file.path(outdir, "example_aoi.gpkg")
aoi_vect <- terra::vect(aoi)
clip_region <- terra::ext(aoi_vect)

# Define hdf5 output file

# Clip the data for only the first ATL08 file
outputs <- atl08_files |> basename() |> map(tempfile, fileext=".h5")
atl08_clippeds <- map2(atl08_h5, outputs, ~clip(.x, output = .y, clip_obj = clip_region))
atl08_clipped_dt <- atl08_clippeds |> map(ATL08_seg_attributes_dt) |> rbindlist2()


head(atl08_clipped_dt)

# Display location of clipped data
atl08_seg_vect <- to_vect(atl08_seg_dt)
atl08_seg_clip_vect <- to_vect(atl08_clipped_dt)

bbox <- terra::vect(terra::ext(atl08_seg_clip_vect), crs = "epsg:4326")
centroid <- terra::geom(terra::centroids(bbox))

map1 <- mapview::mapview(
  atl08_seg_clip_vect,
  col.regions = "yellow",
  alpha.regions = 1,
  lwd = 5,
  map.types = c("Esri.WorldImagery"),
  alpha = 0,
  cex = 2,
  legend = FALSE
)


# Final map
final_map <- map1@map |>
  leaflet::addCircleMarkers(data = atl08_seg_vect, radius = 2) |>
  leaflet::addPolygons(
    data = bbox, fillOpacity = 0, weight = 3, color = "white",
    opacity = 1, dashArray = "5, 1, 0"
  ) |>
  leaflet::addLegend(
    position = "topright",
    colors = c("blue", "yellow", "white"),
    labels = c("atl08_segment", "atl08_segment_clipped", "bbox"),
    opacity = 1
  ) |>
  leaflet::setView(centroid[, "x"][[1]], centroid[, "y"][[1]], zoom = 13)

final_map


###########################
## Clip vector using `id`
###########################
aoi <- file.path(outdir, "example_aoi.gpkg")
aoi_vect <- terra::vect(aoi)


atl08_seg_dt_clip <- ATL08_seg_attributes_dt_clipGeometry(
  atl08_clipped_dt, aoi_vect, split_by = "id"
)

atl08_seg_clip_vect <- to_vect(atl08_seg_dt_clip)

colors <- c("#00FF00", "#FF00FF")

map1 <- mapview::mapview(
  atl08_seg_clip_vect,
  alpha = 0,
  col.regions = colors,
  alpha.regions = 1,
  zcol = "poly_id",
  lwd = 5,
  map.types = c("Esri.WorldImagery"),
  cex = 2,
  legend = FALSE
)

# Final map
final_map <- map1@map |>
  leaflet::addCircleMarkers(data = atl08_seg_vect, color = "blue", radius = 2) |>
  leaflet::addPolygons(
    data = aoi_vect, fillOpacity = 0, weight = 3, color = colors,
    opacity = 1, dashArray = "5, 1, 0"
  ) |>
  leaflet::addLegend(
    position = "topright",
    colors = c("blue", "yellow", colors),
    labels = c("atl08_segment", "atl08_segment_clipped", "aoi_1", "aoi_2"),
    opacity = 1
  ) |>
  leaflet::setView(centroid[, "x"][[1]], centroid[, "y"][[1]], zoom = 13)

final_map




# Generic clip() interface (auto-detects type and dispatches)
# clipped <- clip(data_object, clip_obj = aoi_vect)


# ── Extract Reference Ground Tracks ───────────────────────────────────────────

rgt_line_vect    <- rgt_extract(atl03_h5[[2]])
rgt_polygon_vect <- rgt_extract(atl03_h5[[2]], line = FALSE)

atl08_seg_vect$formatted_h <- atl08_seg_vect$h_canopy
atl08_seg_vect$formatted_h <- sprintf("Height: %.2f", atl08_seg_vect$formatted_h)

map1 <- mapview::mapview(
  atl08_seg_vect,
  zcol = "formatted_h",
  color='#f0f',
  map.types = c("Esri.WorldImagery"),
  cex = 2,
  legend = FALSE
)@map |>
  leaflet::addPolylines(data = rgt_line_vect, color='#ff0', opacity=1) |>
  leaflet::addPolygons(data = rgt_polygon_vect, fill=NA, color='#0f0', opacity=1) |>
  leaflet::addLegend(
    position = "topright",
    colors = c("#f0f", "#ff0", "#0f0"),
    labels = c("ATL08 Segments", "RGT Line", "RGT Polygon"),
    opacity = 1
  )

map1

# ── Join ATL03 and ATL08 ──────────────────────────────────────────────────────

atl03_atl08_dts <- lapply(seq_along(atl03_h5), function(ii) {
  ATL03_ATL08_photons_attributes_dt_join(atl03_h5[[ii]], atl08_h5[[ii]])
})

atl03_atl08_dt <- map2(atl03_h5, atl08_h5, ATL03_ATL08_photons_attributes_dt_join) |>
  rbindlist2()

# Plot photon profiles
layout(matrix(c(1, 2), ncol = 1))
plot(atl03_atl08_dt[orbit_number == 3208], y = "h_ph",
  colors = c("gray", "#bd8421", "forestgreen", "green"),
  xlim = c(32250, 35000), beam = "gt2r", cex = 0.7, pch = 16)

plot(atl03_atl08_dt[orbit_number == 3208], y = "ph_h",
  colors = c("gray", "#bd8421", "forestgreen", "green"),
  xlim = c(32250, 35000), beam = "gt2r", cex = 0.7, pch = 16, legend = NA)

# Raster statistics from joined data
h_canopy_rast <- ATL03_ATL08_photons_attributes_dt_gridStat(
  atl03_atl08_dt[ph_h < 50 & ph_h > 0],
  func = list(h_canopy = quantile(ph_h, 0.98), count = .N),
  res  = 0.01
)

plot(h_canopy_rast)

# ── Custom Segment Length ─────────────────────────────────────────────────────

# Create segments of 30 m instead of the default 100 m
atl03_atl08_photons_grouped_dt <- ATL03_ATL08_segment_create(
  atl03_atl08_dt,
  segment_length = 30,
  centroid       = "mean",
  output         = NA,
  overwrite      = FALSE
)

# Compute statistics per custom segment
atl03_atl08_seg_dt <- ATL03_ATL08_compute_seg_attributes_dt_segStat(
  atl03_atl08_photons_grouped_dt,
  list(
    h_canopy_ge0    = quantile(ph_h, 0.98),
    h_canopy_gt0    = quantile(ph_h[ph_h > 0], 0.98),
    n_ground        = sum(classed_pc_flag == 1),
    n_mid_canopy    = sum(classed_pc_flag == 2),
    n_top_canopy    = sum(classed_pc_flag == 3),
    n_canopy_total  = sum(classed_pc_flag >= 2)
  ),
  ph_class = c(1, 2, 3)
)


# ── Machine Learning: Predict and Rasterize ───────────────────────────────────

library(randomForest)

# Sample training data
h_canopy_train <- c(13.9, 3.1, 2.2, 4.6, 21.6, 7.2, 5, 7.7, 0.8, 9.7,
                    11, 11.3, 15.5, 5.1, 10.4, 0.6, 14.6, 13.3, 9.8, 14.7)
agbd_train     <- c(144.8, 27.5, 51.6, 60.5, 232.3, 102.8, 33.1, 91.3, 23,
                    120.1, 125.7, 127.2, 147.4, 48.8, 103.3, 55.9, 181.8, 139.9, 120.1, 162.8)

set.seed(42)
model <- randomForest(data.frame(h_canopy = h_canopy_train), agbd_train)

# Predict for all ATL08 data and write to HDF5
out_h5 <- tempfile(fileext = ".h5")
for (atl08_h5_item in atl08_h5) {
  seg_dt <- ATL08_seg_attributes_dt(atl08_h5_item, attributes = "h_canopy")
  seg_dt[h_canopy > 100, h_canopy := NA_real_]
  seg_dt <- na.omit(seg_dt)
  predicted_h5 <- predict_h5(model, seg_dt, out_h5)
}

# Rasterize predictions
output_raster <- tempfile(fileext = ".tif")
x    <- predicted_h5[["longitude"]][]
y    <- predicted_h5[["latitude"]][]
bbox <- terra::ext(min(x), max(x), min(y), max(y))
rasterize_h5(predicted_h5, output_raster, bbox = bbox, res = 0.005)

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

###############################################################################
# ── Google Earth Engine Upscaling (requires ee_initialize) ───────────────────

library(terra)
library(leaflet)
library(purrr)
library(ICESat2VegR)

setwd("C:/Users/caioh/Documents/test_ICESat2/")
atl08_files <- list.files(pattern='ATL08.*clip.h5$')
atl08_h5 <- atl08_files |> map(ATL08_read)
atl08_seg_dt <- atl08_h5 |>
  map(ATL08_seg_attributes_dt, attribute = "h_canopy") |>
  rbindlist2()

atl08_seg_dt <- atl08_seg_dt[h_canopy < 20000]

atl08_seg_vect <- to_vect(atl08_seg_dt)
centroid <- geom(centroids(vect(ext(atl08_seg_vect))))


redYellowGreen <- function(n) grDevices::hcl.colors(n, "RdYlGn")
map <- mapview::mapView(
  atl08_seg_vect,
  zcol="h_canopy",
  col.regions = redYellowGreen,
  map.types = c("Esri.WorldImagery")
  )@map |>
  leaflet::setView(lng = centroid[, "x"][[1]], lat = centroid[, "y"][[1]], zoom = 12)

map


# Search and load HLS imagery
hls_search <- search_datasets("Harmonized", "Landsat")
hls_id     <- get_catalog_id(hls_search[1]$id)
hls_collection <- ee$ImageCollection(hls_id)
names(hls_collection)

# Filter, cloud-mask, and compute median composite
### Define area of interest (aoi) clip boundaries and time and cloud mask for filtering.

# For filtering we have the Fmask quality flags from HLS User Guide: <https://lpdaac.usgs.gov/documents/1698/HLS_User_Guide_V2.pdf>
#
#   | Bit number | Mask name                  | Bit value | Mask description        |
#   |------------|----------------------------|-----------|--------------------------|
#   | 7–6        | Aerosol level              | 11        | High aerosol             |
#   | 7–6        | Aerosol level              | 10        | Moderate aerosol         |
#   | 7–6        | Aerosol level              | 01        | Low aerosol              |
#   | 7–6        | Aerosol level              | 00        | Climatology aerosol      |
#   | 5          | Water                      | 1         | Yes                      |
#   | 5          | Water                      | 0         | No                       |
#   | 4          | Snow/ice                   | 1         | Yes                      |
#   | 4          | Snow/ice                   | 0         | No                       |
#   | 3          | Cloud shadow               | 1         | Yes                      |
#   | 3          | Cloud shadow               | 0         | No                       |
#   | 2          | Adjacent to cloud/shadow   | 1         | Yes                      |
#   | 2          | Adjacent to cloud/shadow   | 0         | No                       |
#   | 1          | Cloud                      | 1         | Yes                      |
#   | 1          | Cloud                      | 0         | No                       |
#   | 0          | Cirrus                     | NA        | Reserved (not used)      |


 bitMask <- bitwShiftL(1,1) + bitwShiftL(1,2) + bitwShiftL(1,3)
 bbox <- terra::ext(atl08_seg_vect)

 aoi <- ee$Geometry$BBox(
   west = bbox$xmin,
   south = bbox$ymin,
   east = bbox$xmax,
   north = bbox$ymax
 )

 hls <- hls_collection$
   filterDate("2019-04-01", "2019-05-31")$
   filterBounds(aoi)$
   map(function(x) x$updateMask(!(x[["Fmask"]] & bitMask)))$
   median()


 hls_unmasked <- hls_collection$
   filterDate("2019-04-01", "2019-05-31")$
   filterBounds(aoi)$
   median()


 ###########################
 ## Calculate band indexes
 ###########################
 # Rename bands
 hls_unmasked <- hls_unmasked[["B2", "B3", "B4", "B5", "B6", "B7"]]
 names(hls_unmasked) <- c("blue", "green", "red", "nir", "swir1", "swir2")

 hls <- hls[["B2", "B3", "B4", "B5", "B6", "B7"]]
 names(hls) <- c("blue", "green", "red", "nir", "swir1", "swir2")

 # Add evi
 nir <- hls[["nir"]]
 red <- hls[["red"]]
 blue <- hls[["blue"]]

 hls[["evi"]] <- (2.5 * (nir - red)) / (nir + 6 * red - 7.5 * blue + 1)
 print(hls)


 ## Visualize the resulting image
 centroid <- mean(bbox)
 map <- leaflet::leaflet() |>
   addEEImage(hls, bands = list("red", "green", "blue"), group = "masked", max = 0.6) |>
   addEEImage(hls_unmasked, bands = list("red", "green", "blue"), group = "unmasked", max = 0.6) |>
   setView(lng = centroid[1], lat = centroid[2], zoom = 13) |>
   addLayersControl(
     baseGroups = c("unmasked", "masked"),
     options = layersControlOptions(collapsed = FALSE)
   )

 map


 # Read the map data into segments points

 extracted_dt <- seg_ancillary_extract(hls, atl08_seg_vect)
 extracted_dt


 # Use R to fit randomForest
 library(randomForest)
 bandNames <- names(hls)
 x <- extracted_dt[, .SD, .SDcols = bandNames]
 y <- extracted_dt[["h_canopy"]]

 # Mask the NA values
 na_mask <- y < 100

 x <- x[na_mask]
 y <- y[na_mask]

 set.seed(42)

 rf_model <- randomForest::randomForest(x, y, ntree = 300, mtry = 1)
 print(rf_model)
 rf_importance <- importance(rf_model)
 barplot(rf_importance[, "IncNodePurity"])

 ####
 # Transform to R randomForest to GEE
 gee_model <- build_ee_forest(rf_model)

 # Classify the image
 result <- hls$classify(gee_model)

 ###################
 # View the map

 # Enhance contrast calculating max and min
 min_hcanopy <- min(atl08_seg_dt$h_canopy)
 max_hcanopy <- 20
 atl08_seg_vect$h_canopy <- round(atl08_seg_vect$h_canopy, 3) # Round off to 3 decimal places

 # Color palette
 forest_height_palette <- c("#ffffff", "#99cc99", "#006600", "#004d00")
 atl08_seg_dt$h_canopy <- round(atl08_seg_vect$h_canopy, 3)
 palette_colors <- colorNumeric(forest_height_palette, range(atl08_seg_dt$h_canopy))(unique(atl08_seg_dt[order(h_canopy), h_canopy]))


 modelled_map <- terra::plet(
   atl08_seg_vect,
   "h_canopy",
   col = palette_colors,
   tiles = NULL,
   cex = 2
 ) |>
   addEEImage(
     hls,
     bands = c("red", "green", "blue"),
     group = "hls",
     max = 0.6
   ) |>
   leaflet::addLegend(
     pal = colorNumeric(forest_height_palette, seq(min_hcanopy, max_hcanopy)),
     values = seq(min_hcanopy, max_hcanopy, length = 3),
     opacity = 1,
     title = "h_canopy",
     position = "bottomleft",
   ) |>
   setView(lng = centroid[1], lat = centroid[2], zoom = 12)

 modelled_map


 final <- modelled_map |>
   addEEImage(
     result,
     bands = "classification",
     group = "classification",
     min = min_hcanopy,
     max = max_hcanopy,
     palette = forest_height_palette
   ) |>
   addLayersControl(
     overlayGroups = c("classification"),
     options = layersControlOptions(collapsed = FALSE)
   )

 final
# ── Close Files ───────────────────────────────────────────────────────────────

atl03_h5 |> map(close)
atl08_h5 |> map(close)
