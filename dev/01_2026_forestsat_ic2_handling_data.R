# ICESat2VegR Workflow Script
# Authors: Carlos Alberto Silva and Caio Hamamura
# Processes NASA ICESat-2 ATL03 and ATL08 data for land and vegetation applications

# ── Installation ──────────────────────────────────────────────────────────────
# install.packages('pak')
pak::pkg_install("carlos-alberto-silva/ICESat2VegR", dependencies = TRUE)

# -- Load
library(ICESat2VegR)
library(data.table)

# ── Authentication ────────────────────────────────────────────────────────────
# This is necessary for downloading data from NASA Earthdata
# and for using Google Earth Engine (if needed)

# Option 1: Set env vars EARTHDATA_USERNAME and EARTHDATA_PASSWORD before opening R
# Option 2:
earthdata_login(output_dir = "~")

# Google Earth Engine (if needed)
# Create a project with Earth Engine
# <https://console.cloud.google.com/earth-engine/welcome>
ee_initialize("project-name")


# ── Download Example Dataset ──────────────────────────────────────────────────
outdir <- "~/ic2_workshop_data"
dir.create(outdir, showWarnings = FALSE)

url <- "https://github.com/carlos-alberto-silva/ICESat2VegR/releases/download/example_datasets/Study_Site.zip"
zip_file <- file.path(outdir, "Study_Site.zip")
if (!file.exists(zip_file)) download.file(url, zip_file, mode = "wb")

if (!file.exists(file.path(outdir, "ATL03_20190413122158_02330302_006_02_clip.h5"))) {
  unzip(zip_file, exdir = outdir)
}


# ── Search Parameters ─────────────────────────────────────────────────────────
setwd(outdir)
vector_data <- terra::vect("example_aoi.gpkg")
bounding_box <- terra::ext(vector_data)
bounding_box

lower_left_lon <- bounding_box$xmin
lower_left_lat <- bounding_box$ymin
upper_right_lon <- bounding_box$xmax
upper_right_lat <- bounding_box$ymax
daterange <- c("2019-04-12", "2019-05-03")


# ── Find & Download Granules (Local Workflow) ─────────────────────────────────
atl03_granules_local <- ATLAS_dataFinder(
  short_name = "ATL03",
  lower_left_lon, lower_left_lat, upper_right_lon, upper_right_lat,
  version = "007",
  daterange = daterange,
  persist = FALSE,
  cloud_computing = FALSE
)
atl03_granules_local

atl08_granules_local <- ATLAS_dataFinder(
  short_name = "ATL08",
  lower_left_lon, lower_left_lat, upper_right_lon, upper_right_lat,
  version = "007",
  daterange = daterange,
  persist = TRUE,
  cloud_computing = FALSE
)
atl08_granules_local


# ATLAS_dataDownload(atl03_granules_local, outdir)
# ATLAS_dataDownload(atl08_granules_local, outdir)


# ── Open HDF5 Files ───────────────────────────────────────────────────────────

atl03_files <- list.files(outdir, "ATL03.*h5", full.names = TRUE)
atl03_h5 <- lapply(atl03_files, ATL03_read)

atl08_files <- list.files(outdir, "ATL08.*h5", full.names = TRUE)
atl08_h5 <- lapply(atl08_files, ATL08_read)

# Alternative pipe operator
atl03_h5 |> lapply(close)
atl08_h5 |> lapply(close)

# Alternative purrr
library(purrr)
atl03_h5 <- atl03_files |> map(ATL03_read)
atl08_h5 <- atl08_files |> map(ATL08_read)


# ── Extract ATL03 Photon Attributes ───────────────────────────────────────────
atl03_photons_dt <- atl03_h5 |> map(ATL03_photons_attributes_dt)

atl03_photons_dt[[1]][, name := basename(atl03_files[[1]])]
atl03_photons_dt[[2]][, name := basename(atl03_files[[2]])]

atl03_photons_dt

# Quick plot of a subset
atl03_photons_dt[[2]][beam == "gt1r" & dist_ph_along < 10000] |> with(
  plot(dist_ph_along, h_ph, xlab = "dist_ph_along", ylab = "Elevation (m)", pch = 16, cex = 0.2)
)

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
  ) |>
  rbindlist2()

atl03_seg_dt <- atl03_seg_dt[h_ph < 20000] # remove outliers


# ── Extract ATL08 Segment Attributes ─────────────────────────────────────────

atl08_seg_dt <- atl08_h5 |>
  map(
    ATL08_seg_attributes_dt,
    attributes = c("h_canopy", "h_te_mean", "terrain_slope", "canopy_openness", "night_flag")
  ) |>
  rbindlist2()

atl08_seg_dt <- atl08_seg_dt[h_canopy < 100 & h_te_mean < 20000]

# Histograms
layout(t(1:2))
hist(atl03_seg_dt$h_ph, col = "#bd8421", xlab = "Elevation (m)", main = "ATL03 h_ph")
hist(atl08_seg_dt$h_canopy, col = "green", xlab = "Height (m)", main = "ATL08 h_canopy")
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
  func = max(h_canopy), res = 0.01
)


plot(max_h_canopy, main = "max h_canopy")

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
outputs <- atl08_files |>
  basename() |>
  map(tempfile, fileext = ".h5")
atl08_clippeds <- map2(atl08_h5, outputs, ~ clip(.x, output = .y, clip_obj = clip_region))
atl08_clipped_dt <- atl08_clippeds |>
  map(ATL08_seg_attributes_dt) |>
  rbindlist2()


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
map1

# ── Final Map ──────────────────────────────────────────────────────────────────
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


# ── Clip vector by `id` ───────────────────────────────────────────────────────
aoi <- file.path(outdir, "example_aoi.gpkg")
aoi_vect <- terra::vect(aoi)

atl08_seg_dt_clip <- ATL08_seg_attributes_dt_clipGeometry(
  atl08_clipped_dt, aoi_vect,
  split_by = "id"
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
map1
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

rgt_line_vect <- rgt_extract(atl03_h5[[2]])
rgt_polygon_vect <- rgt_extract(atl03_h5[[2]], line = FALSE)

atl08_seg_vect$formatted_h <- atl08_seg_vect$h_canopy
atl08_seg_vect$formatted_h <- sprintf("Height: %.2f", atl08_seg_vect$formatted_h)

map1 <- mapview::mapview(
  atl08_seg_vect,
  zcol = "formatted_h",
  color = "#f0f",
  map.types = c("Esri.WorldImagery"),
  cex = 2,
  legend = FALSE
)@map |>
  leaflet::addPolylines(data = rgt_line_vect, color = "#ff0", opacity = 1) |>
  leaflet::addPolygons(data = rgt_polygon_vect, fill = NA, color = "#0f0", opacity = 1) |>
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
plot(atl03_atl08_dt[orbit_number == 3208],
  y = "h_ph",
  colors = c("gray", "#bd8421", "forestgreen", "green"),
  xlim = c(32250, 35000), beam = "gt2r", cex = 0.7, pch = 16
)

plot(atl03_atl08_dt[orbit_number == 3208],
  y = "ph_h",
  colors = c("gray", "#bd8421", "forestgreen", "green"),
  xlim = c(32250, 35000), beam = "gt2r", cex = 0.7, pch = 16, legend = NA
)


# Compare with original ATL03 photons
layout(matrix(c(1, 2), ncol = 1))
plot(atl03_atl08_dt[orbit_number == 3208],
  y = "h_ph",
  colors = c("gray", "#bd8421", "forestgreen", "green"),
  xlim = c(32250, 35000), beam = "gt2r", cex = 0.7, pch = 16
)

atl03_photons_dt[
  beam == "gt2r" &
    dist_ph_along < 35000 &
    dist_ph_along > 32250 &
    name == "ATL03_20190413122158_02330302_007_01_clip.h5"
] |> with(
  plot(
    dist_ph_along,
    y = h_ph,
    ylim = c(35, 80),
    xlim = c(32250, 35000), cex = 0.7, pch = 16
  )
)

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
head(atl03_atl08_seg_dt)

# Convert the segments to vector
atl03_atl08_vect <- to_vect(atl03_atl08_seg_dt[h_canopy_gt0 <= 31])

library(leaflet)

centroid <- atl03_atl08_dt[, .(x = mean(lon_ph), y = mean(lat_ph))]
map_output <- mapview::mapView(
  atl03_atl08_vect,
  zcol = "h_canopy_gt0",
  col.regions = grDevices::hcl.colors(9, "RdYlGn"),
  alpha = 0,
  layer.name = "h_canopy",
  map.types = c("Esri.WorldImagery"),
  cex = 4
)@map |>
  leaflet::setView(lng = centroid$x, lat = centroid$y, zoom = 13)

map_output


# ── Close Files ───────────────────────────────────────────────────────────────

atl03_h5 |> map(close)
atl08_h5 |> map(close)
