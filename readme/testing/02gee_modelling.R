library(ICESat2VegR)

geom <- terra::vect("../inst/extdata/all_boundary.shp")
bbox <- terra::ext(geom)
year <- 2019
aprilPlaceholder <- "%s-04-01"
mayPlaceholder <- "%s-05-31"



aoi <- ee$Geometry$BBox(
  west = bbox$xmin,
  south = bbox$ymin,
  east = bbox$xmax,
  north = bbox$ymax
)

aoi <- aoi$buffer(30)

search <- search_datasets("hls")

# get_catalog_id <- get_catalog_path
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
  filterDate(gettextf(aprilPlaceholder, year), gettextf(mayPlaceholder, year))$
  filter("CLOUD_COVERAGE < 10")$
  map(hlsMask)$
  map(waterMask)$
  median()$
  clip(aoi)



hls <- hls[["B2", "B3", "B4", "B5", "B6", "B7"]]
names(hls) <- c("blue", "green", "red", "nir", "swir1", "swir2")

# ndvi
hls[["ndvi"]] <- (hls[["nir"]] - hls[["red"]]) / (hls[["nir"]] + hls[["red"]])



# kndvi
sigma <- 1
nir <- hls[["nir"]]
red <- hls[["red"]]
nir_red <- nir - red
knr <- exp((nir_red)^2 / (2 * sigma^2))
hls[["kndvi"]] <- (1 - knr) / (1 + knr)

# evi
blue <- hls[["blue"]]
hls[["evi"]] <- 2.5 * (nir_red) / (nir + 6 * red - 7.5 * blue + 1)


# savi
hls[["savi"]] <- 1.5 * (nir_red) / (nir + red + 0.5)

# MSAVI2
p1 <- 2 * nir + 1
hls[["msavi"]] <- p1 - sqrt((p1^2) - 8 * (nir_red)) / 2


# Linear spectral unmixing
soil <- c(0.14, 0.16, 0.22, 0.39, 0.45, 0.27)
veg <- c(0.086, 0.062, 0.043, 0.247, 0.109, 0.039)
water <- c(0.07, 0.039, 0.023, 0.031, 0.011, 0.007)

img <- hls[[c("blue", "green", "red", "nir", "swir1", "swir2")]]
unmixing <- img$unmix(list(soil, veg, water))

names(unmixing) <- c("f_soil", "f_veg", "f_water")


hls[[""]] <- unmixing
# Or
# combination <- c(hls, unmixing)

# Simple Ratio Index
hls[["sri"]] <- nir / red

# NDWI
green <- hls[["green"]]
hls[["ndwi"]] <- (green - nir) / (green + nir)

# Green Chloropyl Index
hls[["gci"]] <- (nir / green) - 1

# WDRVI
hls[["wdrvi"]] <- ((0.1 * nir) - red) / ((0.1 * nir) + red)

# Global Vegetation Moisture Index
swir1 <- hls[["swir1"]]
hls[["gvmi"]] <- ((nir + 0.1) - (swir1 + 0.02)) / ((nir + 0.1) + (swir1 + 0.02))

# Chlorophyll Vegetation Index
hls[["cvi"]] <- nir * (red / (green^2))

# Clay minerals ratio
swir2 <- hls[["swir2"]]
hls[["cmr"]] <- swir1 / swir2


# Radar vegetation index
# hls[["rvi"]] <- sqrt(vv/(vv+vh))*(vv/vh)

###########################
# Elevation
###########################
result <- search_datasets("nasa", "dem")
catalog_id <- get_catalog_id(result)

elevation <- ee$Image(catalog_id)

the_slope <- as.integer(slope(as.integer(elevation)) * 1000)
the_aspect <- aspect(elevation)

# stackDem
stackDem <- c(elevation, the_slope, the_aspect)$clip(aoi)


#############################
## SENTINEL 1C
#############################

s1c <- ee$ImageCollection("COPERNICUS/S1_GRD")$
  filterBounds(aoi)$
  filterDate(gettextf(aprilPlaceholder, year), gettextf(mayPlaceholder, year))$
  filter(ee$Filter$listContains("transmitterReceiverPolarisation", "VV"))$
  filter(ee$Filter$listContains("transmitterReceiverPolarisation", "VH"))$
  filter(ee$Filter$eq("instrumentMode", "IW"))

s1c <- s1c$sort("system:time_start", FALSE)
s1c <- s1c$reduce(ee$Reducer$firstNonNull())
s1c <- s1c[["VV_first", "VH_first"]]
names(s1c) <- c("vv", "vh")
s1c <- as.integer(s1c * 100)

# Radar vegetation index
vv <- s1c[["vv"]]
vh <- s1c[["vh"]]

# RVI
s1c[["rvi"]] <- as.integer(sqrt(vv / (vv + vh)) * (vv / vh) * 10000)

# COPOLs
s1c[["copol"]] <- as.integer((vv / vh) * 10000)
s1c[["copol2"]] <- as.integer(((vv - vh) / (vv + vh)) * 10000)
s1c[["copol3"]] <- as.integer((vh / vv) * 10000)



##################
## FINAL STACK
##################
fullStack <- c(hls, s1c, stackDem)

kernel <- ee$Kernel$fixed(3, 3, list(c(1, 1, 1), c(1, 1, 1), c(1, 1, 1)), 3, 3, FALSE)
reducerNames <- c("mean", "min", "max", "stdDev")
for (reducerName in reducerNames) {
  reducer <- ee$Reducer[[reducerName]]()
  fullStack[[""]] <- hls$reduceNeighborhood(reducer, kernel)
  fullStack[[""]] <- stackDem$reduceNeighborhood(reducer, kernel)
}

# Texture
img <- as.integer((hls[[c("blue", "green", "red", "nir", "swir1", "swir2")]]) * 1e4)
fullStack[[""]] <- img$glcmTexture(size = 3)



