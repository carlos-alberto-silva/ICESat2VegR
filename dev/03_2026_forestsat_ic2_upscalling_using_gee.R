###############################################################################
# ── Google Earth Engine Upscaling (requires ee_initialize) ───────────────────

library(terra)
library(leaflet)
library(purrr)
library(ICESat2VegR)

setwd("~/ic2_workshop_data/")
ee_initialize("ee-caiohamamura")
atl08_files <- list.files(pattern = "ATL08.*clip.h5$")
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
  zcol = "h_canopy",
  col.regions = redYellowGreen,
  map.types = c("Esri.WorldImagery")
)@map |>
  leaflet::setView(lng = centroid[, "x"][[1]], lat = centroid[, "y"][[1]], zoom = 12)

map


# Search and load HLS imagery
hls_search <- search_datasets("Harmonized", "Landsat")
hls_id <- get_catalog_id(hls_search[1]$id)
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


bitMask <- bitwShiftL(1, 1) + bitwShiftL(1, 2) + bitwShiftL(1, 3)
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


# # Download (Drive)
out <- map_download(
  ee_image = result, method = "drive",
  region = aoi, scale = 10,
  file_name_prefix = "prediction_2018",
  dsn = file.path(tempdir(), "prediction_2018.tif"),
  drive_folder = "EE_Exports", monitor = TRUE
)
