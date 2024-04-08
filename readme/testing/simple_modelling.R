# devtools::load_all()
# devtools::install_deps()
# devtools::document()
# devtools::install()


library(ICESat2VegR)
library(terra)
library(magrittr)

# ICESat2VegR_configure()

geom <- terra::vect("../inst/extdata/all_boundary.shp")
bbox <- terra::ext(geom)
project_data_dir <- "../inst/extdata"
year <- 2019

##################################
## RETRIEVE ICESAT-2 DATA
##################################
years <- c(2019)
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
  4,
  "FORK"
)
parallel::clusterCall(my.cluster, function() {
  library(ICESat2VegR)
})
parallel::parLapply(my.cluster, all_granules_atl03, function(granule) {
  ATLAS_dataDownload(granule, project_data_dir)
})
parallel::stopCluster(my.cluster)

granules03 <- normalizePath(file.path(project_data_dir, basename(all_granules_atl03)))
granules08 <- normalizePath(file.path(project_data_dir, basename(all_granules_atl08)))
granules03 <- intersect(granules03, gsub("ATL08", "ATL03", granules08))
granules08 <- intersect(granules08, gsub("ATL03", "ATL08", granules03))

if (!all(granules03 == gsub("ATL08", "ATL03", granules08))) {
  stop("Wait, something not right!")
}

target_attributes <- c("h_canopy")
source("readme/testing/00_get_segments_100m.R")


# parallel::stopCluster(my.cluster)

all_data_dt <- data.table::rbindlist(all_data)
nrow(all_data_dt)
prepend_class(all_data_dt, "icesat2.atl08_dt")


# saveRDS(all_data_dt, file = file.path(project_data_dir, "2022atl03_atl08_seg_30m.rds"))
saveRDS(all_data_dt, file = file.path(project_data_dir, gettextf("%datl08_seg_20m_100m.rds", year)))
#all_data_dt <- readRDS(file.path(project_data_dir, "2019atl03_atl08_seg_30m.rds"))

################################
## GET EARTH ENGINE STACK
################################
source("readme/testing/02gee_modelling.R")

###########################
## MODELLING
###########################
set.seed(47289143)
degree_to_meter_factor <- 111139
out_dt2 <- out_dt[night_flag == 1]
atl08_20m <- rbindlist(list(
  out_dt2[, list(
    h_canopy = h_canopy_20m_geo_1,
    longitude = longitude_20m_1,
    latitude =  latitude_20m_1
  )],
  out_dt2[, list(
    h_canopy = h_canopy_20m_geo_2,
    longitude = longitude_20m_2,
    latitude =  latitude_20m_2
  )],
  out_dt2[, list(
    h_canopy = h_canopy_20m_geo_3,
    longitude = longitude_20m_3,
    latitude =  latitude_20m_3
  )],
  out_dt2[, list(
    h_canopy = h_canopy_20m_geo_4,
    longitude = longitude_20m_4,
    latitude =  latitude_20m_4
  )],
  out_dt2[, list(
    h_canopy = h_canopy_20m_geo_5,
    longitude = longitude_20m_5,
    latitude =  latitude_20m_5
  )]
))

atl08_20m_filter <- atl08_20m[h_canopy < 50]
prepend_class(atl08_20m_filter, "icesat2.atl08_dt")

sampled_vect <- to_vect(atl08_20m_filter)
sampled_vect$I <- seq_along(sampled_vect)
out_dt <- seg_gee_ancillary_dt_extract(fullStack, sampled_vect, scale = 30, chunk_size = 2000)
# out_dt2 <- seg_gee_ancillary_dt_extract(fullStack, sampled_vect[22001:24000], scale = 30, chunk_size = 2000)



cols_add <- setdiff(colnames(all_data_dt), colnames(out_dt))
out_dt[, I(cols_add) := all_data_dt[out_dt$I, .SD, .SDcols = cols_add]]
prepend_class(out_dt, "icesat2.atl03_atl08_seg_dt")
class(out_dt)
saveRDS(out_dt, file.path(project_data_dir, gettextf("%dextracted_dt_20m.rds", year)))


out_dt2 <- out_dt[n_canopy_total > 5 & n_ground > 5]
{
sampled <- sample(all_data_dt2, method = spacedSampling(0.99999999, radius = 30 / degree_to_meter_factor, spatialIndex  = index))
sampled2 <- sample(all_data_dt2[sampled[, I]], stratifiedSampling(0.7, "h_canopy_ge0", breaks = c(0, 10, 20, 999)))

library(ggplot2)
png("here.png")
ggplot() + 
  geom_point(data = sampled2, aes(longitude, latitude, color = h_canopy_ge0)) +
  scale_color_gradientn(colours=c("red", "yellow", "forestgreen", "forestgreen", "forestgreen", "forestgreen", "forestgreen", "forestgreen"), breaks = c(0, 15, 30, 99), limits = c(0, 99))
}
dev.off()
set.seed(47289143)

out_dt2 = out_dt[night_flag == 1]
xcols <- names(out_dt)[-c(1:3, (ncol(out_dt) - 21):ncol(out_dt))]


y <- "h_canopy"

nboots <- 10

res <- var_select_local(out_dt, y, xcols = xcols, "forward", nboots = nboots)

print(res)
selected_properties <- res$property[nrow(res)-1]
selected_cols <- strsplit(selected_properties, ", ")[[1]]
out_dt2[, I := .I]
train_x <- sample(out_dt2, randomSampling(0.7))
test <- out_dt2[-train_x$I]

rf_model <- randomForest::randomForest(train_x[, .SD, .SDcols = selected_cols], train_x[[y]], method = "randomForest")
class(rf_model)

predicted <- predict(rf_model, test)
png("here.png")
stats_model(test[[y]], predicted, xlim = c(0, 35), ylim = c(0, 35))
dev.off()

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
