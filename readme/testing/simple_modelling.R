# devtools::install_deps()
# devtools::install()

library(ICESat2VegR)
library(terra)
library(magrittr)

geom <- terra::vect("Z:/01_Projects/04_NASA_ICESat2/00_data/01_SHPs/all_boundary.shp")
bbox <- terra::ext(geom)

##################################
## RETRIEVE ICESAT-2 DATA
##################################
years <- c(2019:2022)
aprilPlaceholder <- "%s-04-01"
mayPlaceholder <- "%s-05-31"

all_granules <- c()

for (year in years) {
    granules <- ICESat2VegR::ICESat2_finder(
        short_name = "ATL08",
        version = "006",
        daterange = c(gettextf(aprilPlaceholder, year), gettextf(mayPlaceholder, year)),
        lower_left_lon = bbox$xmin,
        lower_left_lat = bbox$ymin,
        upper_right_lon = bbox$xmax,
        upper_right_lat = bbox$ymax
    )
    all_granules <- c(all_granules, granules)
}

ICESat2_download(all_granules, "Z:\\01_Projects\\04_NASA_ICESat2\\00_data\\04_ICESat2_datasets\\StudySite")
granules <- list.files("Z:\\01_Projects\\04_NASA_ICESat2\\00_data\\04_ICESat2_datasets\\StudySite", "ATL08.*h5", full.names = TRUE)

power_beams <- c(
    "gt1r",
    "gt2r",
    "gt3r"
)

target_attributes <- c("h_canopy")

all_dt <- list()
ii <- 1
for (granule in granules) {
    atl08 <- ATL08_read(granule)
    all_dt[[ii]] <- ATL08_seg_attributes_dt(atl08, beam = power_beams, attribute = target_attributes)[h_canopy < 100]
    ii <- ii + 1
}
dt <- data.table::rbindlist(all_dt)
dt
prepend_class(dt, "icesat2.atl08_dt")

dt2 <- ATL08_seg_attributes_dt_clipGeometry(dt, polygon = geom)
nrow(dt2)


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

search <- search_datasets("hls", "landsat")
id <- get_catalog_id(search)
collection <- ee$ImageCollection(id)


cloudMask <-
    2^0 +
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
    filterDate("2020-02-01", "2020-05-31")$
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
sampled <- sample(dt2, method = spacedSampling(100, radius = 30 / degree_to_meter_factor))

sampled_vect <- terra::vect(
    as.data.frame(sampled[, .SD, .SDcols = c("beam", "longitude", "latitude", "h_canopy")]),
    geom = c("longitude", "latitude")
)

out_dt <- seg_gee_ancillary_dt_extract(fullStack, sampled_vect, scale = 30)

x <- out_dt[, .SD, .SDcols = names(fullStack)]
y <- out_dt[, "h_canopy"]
nboots <- 3

res <- var_select(x, y, "forward", nboots = nboots)
print(res)
selected_properties <- res$properties[-nrow(res)]
x_selected <- x[, .SD, .SDcols = selected_properties]

rf_model <- model_fit(x_selected, y, method = "randomForest")
predicted <- rf_model %>% predict(x)

stat_model(y[[1]], predicted)