# ======================================================================
# Upscaling ICESat-2 canopy height using AlphaEarth Embeddings and ancillary datasets
# ======================================================================


# ======================================================================
# Introduction
# ======================================================================


# ======================================================================
# Install and load required packages
# ======================================================================

# The r-universe version (recommended for the latest version)
install.packages('ICESat2VegR', repos = c('https://carlos-alberto-silva.r-universe.dev', 'https://cloud.r-project.org'))
# The CRAN version
install.packages("ICESat2VegR")

# Required additional libraries
need_pkgs <- c(                                                                       # Packages required by the workflow
  "reticulate",   # Python <-> R interface
  "leaflet",      # interactive maps
  "sf",           # spatial vectors
  "terra",        # rasters & vectors
  "data.table",   # fast tables
  "dplyr"         # tidy helpers
)
missing <- need_pkgs[!need_pkgs %in% rownames(installed.packages())]                 # Identify missing packages
if (length(missing)) {                                                                # If any are missing
  message("Installing missing R packages: ", paste(missing, collapse = ", "))         # Inform the user
  install.packages(missing, repos = repos, dependencies = TRUE)                       # Install from repos with dependencies
}
if (!requireNamespace("mapview", quietly = TRUE)) {
  install.packages("mapview", repos = "https://cloud.r-project.org")
}

library(mapview)
suppressPackageStartupMessages({                                                       # Suppress startup messages for clean logs
  library(ICESat2VegR)                                                                 # ICESat-2 vegetation tools
  library(reticulate)                                                                  # Interface to Python
  library(leaflet)                                                                     # Interactive maps
  library(sf)                                                                          # Simple Features for vector data
  library(terra)                                                                       # Raster + vector geospatial ops
  library(data.table)                                                                  # Fast data tables
  library(dplyr)
  # Tidy verbs
})
pkgs <- c("mapview", "caret")

for (pkg in pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
  library(pkg, character.only = TRUE)
}
cat("Loading libraries... done.\n")                                                    # Confirm loading



# ======================================================================
# Read AOI and define date range
# ======================================================================

aoi_path <- system.file("extdata", "aoi_4326.geojson", package = "ICESat2VegR")
boundary <- sf::st_read(aoi_path, quiet = TRUE)
box      <- sf::st_bbox(boundary)

lower_left_lon  <- box["xmin"];  lower_left_lat  <- box["ymin"]
upper_right_lon <- box["xmax"];  upper_right_lat <- box["ymax"]

daterange <- c("2025-08-01", "2025-08-31")


# ======================================================================
# Discover ATL03 and ATL08 granules
# ======================================================================

atl03_granules_cloud <- ICESat2VegR::ATLAS_dataFinder(
  short_name = "ATL03",
  lower_left_lon, lower_left_lat, upper_right_lon, upper_right_lat,
  version = "007", daterange = daterange, persist = TRUE, cloud_computing = FALSE
)

atl08_granules_cloud <- ICESat2VegR::ATLAS_dataFinder(
  short_name = "ATL08",
  lower_left_lon, lower_left_lat, upper_right_lon, upper_right_lat,
  version = "007", daterange = daterange, persist = TRUE, cloud_computing = FALSE
)

# ======================================================================
# Download granules locally
# ======================================================================

outdir <- tempdir()
ATLAS_dataDownload(atl03_granules_cloud, outdir)
ATLAS_dataDownload(atl08_granules_cloud, outdir)

# ======================================================================
# Pair ATL03 and ATL08 files by shared timestamp
# ======================================================================

atl03_files <- list.files(outdir, pattern = "ATL03.*h5", full.names = TRUE)
atl08_files <- list.files(outdir, pattern = "ATL08.*h5", full.names = TRUE)
get_timestamp <- function(f) sub(".*_(\\d{14})_.*", "\\1", basename(f))
timestamps_comuni <- intersect(sapply(atl08_files, get_timestamp),
                               sapply(atl03_files, get_timestamp))

# ======================================================================
# Read one ATL03/ATL08 pair and join photon attributes
# ======================================================================

ts    <- timestamps_comuni[1]
atl03 <- ICESat2VegR::ATL03_read(atl03_files[grepl(ts, atl03_files)])
atl08 <- ICESat2VegR::ATL08_read(atl08_files[grepl(ts, atl08_files)])
dt <- ICESat2VegR::ATL03_ATL08_photons_attributes_dt_join(atl03, atl08)
head(dt)

# ======================================================================
# Create 20 m segments and compute canopy metrics
# ======================================================================

seg20 <- ATL03_ATL08_segment_create(dt, 20, centroid = "mean", output = NA, overwrite = FALSE)
stats20 <- ATL03_ATL08_compute_seg_attributes_dt_segStat(
  seg20,
  list(
    rh98           = quantile(ph_h[ph_h > 0 & classed_pc_flag %in% c(2,3)], 0.98),
    n_ground       = sum(classed_pc_flag == 1),
    n_top_canopy   = sum(classed_pc_flag == 3),
    n_canopy_total = sum(classed_pc_flag >= 2),
    mean_solar     = mean(solar_elevation),
    night_flag2    = as.integer(mean(night_flag) > 0.5)
  ),
  ph_class = c(2, 3)
)
# Convert to SpatVector first
stats20_vect <- to_vect(stats20)
centroid <- stats20[, .(x = mean(longitude), y = mean(latitude))]
map_output <- mapview::mapview(
  stats20_vect,
  zcol        = "rh98",
  col.regions = grDevices::hcl.colors(9, "RdYlGn"),
  alpha       = 0,
  layer.name  = "rh98 (m)",
  map.types   = c("Esri.WorldImagery"),
  cex         = 4
)@map %>%
  leaflet::setView(lng = centroid$x, lat = centroid$y, zoom = 10)

map_output


# ======================================================================
# Clip segments to AOI and export GeoJSON
# ======================================================================

stats20_clip <- ATL03_ATL08_seg_attributes_dt_clipGeometry(stats20, boundary)
stats20_clip$year <- as.integer(substr(ts, 1, 4))
stats20_clip <- stats20_clip[stats20_clip$rh98 <= 50]

stats20_clip_sf <- st_as_sf(
  as.data.frame(stats20_clip),
  coords = c("longitude", "latitude"), crs = 4326, remove = FALSE
)
st_write(stats20_clip_sf,
         file.path(outdir, paste0("ATL03_ATL08_", ts, "_20.geojson")),
         delete_dsn = TRUE, quiet = TRUE)
# ======================================================================
# Package configuration
# ======================================================================

Sys.setenv(EE_PROJECT = "your-ee-project")
ICESat2VegR::ee_initialize()

boundary_simple <- sf::st_simplify(
  sf::st_union(sf::st_make_valid(boundary)),
  dTolerance = 0.001
)
region <- .as_ee_geom(boundary_simple)


# ======================================================================
# Build AlphaEarth predictor stack and extract values at segment locations
# ======================================================================

seg_path <- system.file("extdata", "ATL03_ATL08_example_segments.geojson", package = "ICESat2VegR")
data_raw <- sf::read_sf(seg_path)
df_dt    <- data.table::as.data.table(sf::st_drop_geometry(data_raw))
class(df_dt) <- c("icesat2.atl03_atl08_seg_dt", "data.table", "data.frame")

set.seed(1)
df_sampled <- ICESat2VegR::sample(df_dt, method = randomSampling(500))

df_vect <- terra::vect(
  as.data.frame(df_sampled), geom = c("longitude", "latitude"), crs = "EPSG:4326"
)

yr          <- 2025
predictors2 <- ee_build_AlphaEarth_embedding_terrain_stack(region, yr, yr)
sampling_df <- ICESat2VegR::seg_ancillary_extract(predictors2, df_vect, 20)
head(sampling_df[, 1:12])


# ======================================================================
# Visualize the predictor stack as false-color RGB
# ======================================================================

predictors2_rgb <- predictors2$select(c("A00", "A20", "A40"))
#
centroid_lon <- mean(terra::ext(terra::vect(boundary))[c(1,2)])
centroid_lat <- mean(terra::ext(terra::vect(boundary))[c(3,4)])
#
leaflet::leaflet() |>
  addEEImage(
    predictors2_rgb,
    bands = c("A00", "A20", "A40"),
    group = "RGB Embedding",
    min = c(-0.06),
    max = c( 0.12)
  ) |>
  leaflet::addControl(
    html     = "<b>RGB Embedding</b><br>Bands: A00 / A20 / A40",
    position = "bottomleft"
  ) |>
  leaflet::setView(lng = centroid_lon, lat = centroid_lat, zoom = 16) |>
  leaflet::addLayersControl(
    overlayGroups = "RGB Embedding",
    options = leaflet::layersControlOptions(collapsed = FALSE)
  )


# ======================================================================
# Variable selection with RFE
# ======================================================================
set.seed(123)
x <- sampling_df %>%
  dplyr::select(starts_with("A"), "slope", "aspect", "elevation") %>%
  data.frame()
y <- sampling_df %>% dplyr::select("rh98") %>% data.frame()

sel_rf_rfe   <- varSel(x, y$rh98)
best_metrics <- sel_rf_rfe$selvars
print(best_metrics)

# Build importance tables from sel_rf_rfe
imp_desc <- sel_rf_rfe$importance %>%
  dplyr::mutate(selected = parameter %in% sel_rf_rfe$selvars) %>%
  dplyr::arrange(importance)

best_imp_rfe <- imp_desc %>%
  dplyr::filter(selected == TRUE) %>%
  dplyr::arrange(importance)

par(mfrow = c(1, 2), mar = c(4, 6, 2, 1))

# Left — all RFE evaluated, colored by selection
barplot(
  rev(imp_desc$importance),
  names.arg = rev(as.character(imp_desc$parameter)),
  horiz     = TRUE,
  las       = 1,
  main      = "RFE importance\n(green = selected)",
  xlab      = "Importance",
  col       = ifelse(rev(imp_desc$selected), "#1B7837", "grey70"),
  cex.names = 0.7
)
abline(v = 0.2, lty = 2, col = "red")

# Right — selected only
barplot(
  best_imp_rfe$importance,
  names.arg = as.character(best_imp_rfe$parameter),
  horiz     = TRUE,
  las       = 1,
  main      = "Selected predictors",
  xlab      = "Importance",
  col       = "#1B7837",
  cex.names = 0.8
)
abline(v = 0.2, lty = 2, col = "red")

par(mfrow = c(1, 1), mar = c(5, 4, 4, 2))

# ======================================================================
# [1] "A07"  "A22"  "A24"  "A40"  "A56"  "A62"  "A36"  "A34"  "A38"
# ======================================================================
# ======================================================================
# [10] "elevation" "A23" "A20"  "A18"  "A02"  "A57"  "A13"  "A21"  "A30"
# ======================================================================
# ======================================================================
# Train/test split and fit Random Forest model
# ======================================================================

x      <- dplyr::select(sampling_df, dplyr::any_of(best_metrics))
y      <- sampling_df$rh98
data_i <- data.frame(y, x, check.names = FALSE)

idx_train <- caret::createDataPartition(y = data_i$y, p = 0.70, list = FALSE)
trainData <- data_i[ idx_train, ]
testData  <- data_i[-idx_train, ]

fit_rf   <- fit_model(x = dplyr::select(trainData, dplyr::any_of(best_metrics)),
                      y = trainData$y,
                      rf_args = list(ntree = 100))
rf_model <- fit_rf$model

pred_train <- predict(rf_model, newdata = dplyr::select(trainData, dplyr::any_of(best_metrics)))
pred_test  <- predict(rf_model, newdata = dplyr::select(testData,  dplyr::any_of(best_metrics)))

results_train <- fit_metrics(trainData$y, as.numeric(pred_train))
results_test  <- fit_metrics(testData$y,  as.numeric(pred_test))

head(results_train)
head(results_test)

# ======================================================================
# Wall-to-wall canopy height map in GEE
# ======================================================================

mosaic <- ee_build_AlphaEarth_embedding_terrain_stack(region, yr, yr)

ch_map <- map_create(
  model    = fit_rf$model,
  stack    = mosaic,
  aoi      = region,
  reducer  = "mosaic",
  mode     = "auto",
  to_float = TRUE
)

# ======================================================================
# Visualize the canopy height map
# ======================================================================

min_hcanopy <- 0
max_hcanopy <- 40

forest_palette <- colorRampPalette(
  c("#fcd88f", "#00441B", "#1B7837", "#A6DBA0", "#E7E1EF", "#762A83")
)(128)

pal_fun <- leaflet::colorNumeric(forest_palette, domain = c(min_hcanopy, max_hcanopy))

leaflet::leaflet() |>
  addEEImage(
    ch_map,
    bands   = "prediction_layer",
    group   = "Canopy Height",
    min     = min_hcanopy,
    max     = max_hcanopy,
    palette = forest_palette
  ) |>
  leaflet::addLegend(
    pal      = pal_fun,
    values   = c(min_hcanopy, max_hcanopy),
    opacity  = 1,
    title    = "Canopy Height (m)",
    position = "bottomleft"
  ) |>
  leaflet::setView(lng = centroid_lon, lat = centroid_lat, zoom = 15) |>
  leaflet::addLayersControl(
    overlayGroups = "Canopy Height",
    options = leaflet::layersControlOptions(collapsed = FALSE)
  )

# ======================================================================
# Export map to GeoTIFF via Google Drive
# ======================================================================

googledrive::drive_auth(cache = FALSE)

out <- map_download(
  ee_image         = ch_map,
  method           = "drive",
  region           = region,
  scale            = 10,
  file_name_prefix = "prediction_2025",
  dsn              = file.path(outdir, "prediction_2025.tif"),
  drive_folder     = "EE_Exports",
  monitor          = TRUE
)
terra::plot(out)

# ======================================================================
# Close the files
# ======================================================================
