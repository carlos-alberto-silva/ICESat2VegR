# devtools::load_all()
# devtools::install_deps()
# devtools::document()
# devtools::install()


library(ICESat2VegR)
# library(terra)
# library(magrittr)

# ICESat2VegR_configure()

geom <- terra::vect("Z:/01_Projects/04_NASA_ICESat2/00_data/01_SHPs/all_boundary.shp")
bbox <- terra::ext(geom)
project_data_dir <- "Z:\\01_Projects\\04_NASA_ICESat2\\00_data\\04_ICESat2_datasets\\StudySite"
year <- 2022

##################################
## RETRIEVE ICESAT-2 DATA
##################################
years <- c(2022)
aprilPlaceholder <- "%s-04-01"
mayPlaceholder <- "%s-05-31"

all_granules_atl08 <- c()

for (year in years) {
  granules <- ICESat2VegR::ATLAS_dataFinder(
    short_name = "ATL08",
    version = "006",
    daterange = c(gettextf(aprilPlaceholder, year), gettextf(mayPlaceholder, year)),
    lower_left_lon = bbox$xmin,
    lower_left_lat = bbox$ymin,
    upper_right_lon = bbox$xmax,
    upper_right_lat = bbox$ymax
  )
  all_granules_atl08 <- c(all_granules_atl08, granules)
}

ATLAS_dataDownload(all_granules_atl08, project_data_dir)

library(parallel)
library(foreach)
# install.packages("doParallel")
library(doParallel)


all_granules_atl03 <- c()
for (year in years) {
  granules <- ICESat2VegR::ATLAS_dataFinder(
    short_name = "ATL03",
    version = "006",
    daterange = c(gettextf(aprilPlaceholder, year), gettextf(mayPlaceholder, year)),
    lower_left_lon = bbox$xmin,
    lower_left_lat = bbox$ymin,
    upper_right_lon = bbox$xmax,
    upper_right_lat = bbox$ymax
  )
  all_granules_atl03 <- c(all_granules_atl03, granules)
}

my.cluster <- parallel::makeCluster(
  4
)
doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()
foreach::getDoParWorkers()
parallel::clusterCall(my.cluster, function() {
  devtools::load_all()
})

foreach(granule = all_granules_atl03) %dopar% {
  ATLAS_dataDownload(granule, project_data_dir)
}




granules03 <- normalizePath(file.path(project_data_dir, basename(all_granules_atl03)))
granules08 <- normalizePath(file.path(project_data_dir, basename(all_granules_atl08)))
granules03 <- intersect(granules03, gsub("ATL08", "ATL03", granules08))
granules08 <- intersect(granules08, gsub("ATL03", "ATL08", granules03))

if (!all(granules03 == gsub("ATL08", "ATL03", granules08))) {
  stop("Wait, something not right!")
}

target_attributes <- c("h_canopy")
all_data <- list()


# stop(my.cluster)

foreach(ii = seq_along(granules03), .export = c("granules03", "granules08")) %do% {
  atl03_h5 <- ATL03_read(granules03[ii])
  atl08_h5 <- ATL08_read(granules08[ii])
  atl03_atl08_dt <- ATL03_ATL08_photons_attributes_dt_join(atl03_h5, atl08_h5, power_beam_filter = TRUE)
  out_name <- gsub("ATL03", "ATL03_ATL08", basename(granules03[ii]))

  atl03_atl08_dt <- atl03_atl08_dt[
    ph_h < 100 & 
    night_flag == 1 & 
    ph_h >= 0 
    ]
  if (nrow(atl03_atl08_dt) > 0) {
    atl03_atl08_dt_seg <- ATL03_ATL08_segment_create(
      atl03_atl08_dt,
      segment_length = 30
    )

    atl03_atl08_seg_dt <- ATL03_ATL08_compute_seg_attributes_dt_segStat(
      atl03_atl08_dt_seg,
      list(
        h_canopy_ge0 = quantile(ph_h, 0.98),
        h_canopy_gt0 = quantile(ph_h[ph_h > 0], 0.98),
        n_ground = sum(classed_pc_flag == 1),
        n_mid_canopy = sum(classed_pc_flag == 2),
        n_top_canopy = sum(classed_pc_flag == 3),
        n_canopy_total = sum(classed_pc_flag >= 2)
      ),
      ph_class = c(1, 2, 3)
    )

    prepend_class(atl03_atl08_seg_dt, "icesat2.atl08_dt")
    atl03_atl08_seg_dt <- ATL08_seg_attributes_dt_clipGeometry(atl03_atl08_seg_dt, geom)
    atl03_atl08_seg_dt$nid <- NULL

    if (nrow(atl03_atl08_seg_dt) > 0) {
      all_data[[out_name]] <- atl03_atl08_seg_dt
    }
  }

  close(atl03_h5)
  close(atl08_h5)
}

# parallel::stopCluster(my.cluster)
all_data_dt <- data.table::rbindlist(all_data)

# saveRDS(all_data_dt, file = file.path(project_data_dir, "2022atl03_atl08_seg_30m.rds"))

all_data_dt <- readRDS(file.path(project_data_dir, "2022atl03_atl08_seg_30m.rds"))
nrow(all_data_dt)
all_data_dt[h_canopy_ge0 != h_canopy_gt0]
prepend_class(all_data_dt, "icesat2.atl08_dt")

################################
## GET EARTH ENGINE STACK
################################
aoi <- ee$Geometry$BBox(
  west = bbox$xmin,
  south = bbox$ymin,
  east = bbox$xmax,
  north = bbox$ymax
)

aoi <- aoi$buffer(30)

search <- search_datasets("hls")
id <- get_catalog_id(search)
collection <- ee$ImageCollection(id)


cloudMask <- 2^0 +
  2^1 +
  2^2

hlsMask <- function(image) {
  image$updateMask(
    !(image[["Fmask"]] & cloudMask)
  )
}

waterMask <- function(image) {
  image$updateMask(
    image[["B5"]] >= 0.2
  )
}


hls <- collection$
  filterBounds(aoi)$
  filterDate(gettextf(aprilPlaceholder, year), gettextf(mayPlaceholder, year))$
  filter("CLOUD_COVERAGE < 10")$
  map(hlsMask)$
  map(waterMask)$
  median()$
  clip(aoi)



hls <- hls[["B2", "B3", "B4", "B5", "B6", "B7"]]
names(hls) <- c("blue", "green", "red", "nir", "swir1", "swir2")


fullStack <- hls

###########################
## MODELLING
###########################
set.seed(47289143)
degree_to_meter_factor <- 111139

sampled_vect <- to_vect(all_data_dt)
sampled_vect$I <- seq_along(sampled_vect)
out_dt <- seg_gee_ancillary_dt_extract(fullStack, sampled_vect, scale = 30, chunk_size = 1000)

cols_add <- setdiff(colnames(all_data_dt), colnames(out_dt))
out_dt[, I(cols_add) := all_data_dt[out_dt$I, .SD, .SDcols = cols_add]]
prepend_class(out_dt, "icesat2.atl03_atl08_seg_dt")
class(out_dt)
saveRDS(out_dt, file.path(project_data_dir, "2022extracted_dt.rds"))

index <- ANNIndex$new(out_dt$longitude, out_dt$latitude)

nrow(all_data_dt)
nrow(out_dt)


out_dt2 <- out_dt[n_canopy_total > 3 & n_ground > 3]
{
sampled <- sample(out_dt2, method = spacedSampling(0.99999999, radius = 30 / degree_to_meter_factor, spatialIndex  = index))
sampled2 <- sample(all_data_dt[sampled[, I]], stratifiedSampling(0.7, "h_canopy_ge0", breaks = c(0, 10, 20, 999)))

library(ggplot2)

ggplot() + 
  geom_point(data = sampled2, aes(longitude, latitude, color = h_canopy_ge0)) +
  scale_color_gradientn(colours=c("red", "yellow", "forestgreen", "forestgreen", "forestgreen", "forestgreen", "forestgreen", "forestgreen"), breaks = c(0, 15, 30, 99), limits = c(0, 99))
}

set.seed(47289143)
bandNames <- names(fullStack)
x <- out_dt2[, .SD, .SDcols = bandNames]
y <- out_dt2[["h_canopy_ge0"]]

x <- out_dt[, .SD, .SDcols = names(fullStack)]
y <- out_dt[, "h_canopy"]
nboots <- 3

res <- var_select(x, y, "forward", nboots = nboots)
print(res)
selected_properties <- res$properties[-nrow(res)]
x_selected <- x[, .SD, .SDcols = selected_properties]

rf_model <- model_fit(x_selected, y, method = "randomForest")
class(rf_model)

predicted <- predict(rf_model, x)
stats_model(y[[1]], predicted, xlim = c(0, 35), ylim = c(0, 35))


result <- map_create(rf_model, fullStack)


######################
## VISUALIZE THE MAP
######################

forest_height_palette <- c("#ffffff", "#8b4513", "#99cc99", "#006600", "#004d00")


vis <- result$visualize(
  bands = "classification",
  min = 0,
  max = 30,
  palette = forest_height_palette
)
url <- vis$getMapId()$tile_fetcher$url_format

library(leaflet)
center <- apply(matrix(bbox, nrow = 2), 2, mean)
coords <- terra::geom(sampled_vect)
leaflet_map <- leaflet::leaflet() %>%
  addProviderTiles(providers$Esri.WorldImagery, group = "Other") %>%
  leaflet::addTiles(
    urlTemplate = url,
    options = leaflet::tileOptions(opacity = 1),
    group = "h_canopy"
  ) %>%
  leaflet::addCircleMarkers(
    lng = coords[, "x"],
    lat = coords[, "y"],
    radius = 2,
    stroke = FALSE,
    fillOpacity = 1,
    fillColor = "yellow"
  ) %>%
  leaflet::addLegend(
    position = "bottomright",
    pal = colorNumeric(forest_height_palette, 0:30),
    values = seq(0, 30, length = 3),
    opacity = 1,
    title = "H_canopy"
  ) %>%
  addLayersControl(
    overlayGroups = c("h_canopy"),
    options = layersControlOptions(collapsed = FALSE)
  ) %>%
  leaflet::setView(lng = center[1], lat = center[2], zoom = 11)
leaflet_map
