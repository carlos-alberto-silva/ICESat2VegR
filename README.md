---
output: github_document
---

![](https://github.com/carlos-alberto-silva/ICESat2VegR/blob/master/readme/cover.png)<br/>
[![R-CMD-check](https://github.com/carlos-alberto-silva/ICESat2VegR/actions/workflows/r.yml/badge.svg?branch=master)](https://github.com/carlos-alberto-silva/ICESat2VegR/actions/workflows/r.yml)
[![CRAN](https://www.r-pkg.org/badges/version/ICESat2VegR)](https://cran.r-project.org/package=ICESat2VegR)
![Github](https://img.shields.io/badge/Github-0.1.12-green.svg)
![licence](https://img.shields.io/badge/Licence-GPL--3-blue.svg) 
![Downloads](https://cranlogs.r-pkg.org/badges/grand-total/ICESat2VegR)
[![Build Status](https://travis-ci.com/carlos-alberto-silva/ICESat2VegR.svg?token=Jqizwyc6gBxNafNccTdU&branch=master)](https://travis-ci.com/carlos-alberto-silva/ICESat2VegR)

**ICESat2VegR: An R Package for NASA's Ice, Cloud, and Elevation Satellite (ICESat-2) Data Processing and Visualization for Land and Vegetation Applications.**

Authors: Carlos Alberto Silva and Caio Hamamura  

The ICESat2VegR package provides functions for downloading, reading, visualizing, processing and exporting 
NASA's ICESat-2 ATL03 (Global Geolocated Photon Data) and ATL08 (Land and Vegetation Height) 
products for Land and Vegetation Applications in R environment.

# Getting started
















```r
# The CRAN version
install.packages("ICESat2VegR")

# Development version
remotes::install_github("https://github.com/carlos-alberto-silva/ICESat2VegR")
```






## Load the package


```r
library(ICESat2VegR)
```


## Configuring the package

This package uses three Python packages through `reticulate`:

1. [earthaccess](https://github.com/nsidc/earthaccess): allows reading directly from the cloud
2. [h5py](https://github.com/h5py/h5py): for reading hdf5 content from the cloud
3. [earthengine-api](https://github.com/google/earthengine-api): integration with Google Earth Engine for sampling and extracting raster data and upscalling models.

For configuring the package you can use:


```r
ICESat2VegR_configure()

```

This will install miniconda if not available and the necessary packages.

### Notes

 - There are some issues regarding some Python packages not being compatible with the Python version. The above configure function will also try to update python version in that case. 
 - The configure function will warn you about the need to restart R





## Introduction

There are two different ways of working with the ICESat-2 data. Locally or using cloud computing. Common users should work locally, unless they are working within an AWS cloud computing within zone us-west-2.





## Search parameters


```r
# Specifying bounding box coordinates
lower_left_lon <- -96.0
lower_left_lat <- 40.0
upper_right_lon <- -100
upper_right_lat <- 42.0


# Specifying the date range
daterange <- c("2021-10-02", "2021-10-03")
```

## Working locally

First we need to find the granules:


```r
atl03_granules_local <- ATLAS_dataFinder(
  short_name = "ATL03",
  lower_left_lon,
  lower_left_lat,
  upper_right_lon,
  upper_right_lat,
  version = "006",
  daterange = daterange,
  persist = TRUE,
  cloud_computing = FALSE
)

head(atl03_granules_local)
#>      C2596864127-NSIDC_CPRD                                                                                                                      
#> [1,] "https://data.nsidc.earthdatacloud.nasa.gov/nsidc-cumulus-prod-protected/ATLAS/ATL03/006/2021/10/02/ATL03_20211002001658_01461302_006_01.h5"
#> [2,] "https://data.nsidc.earthdatacloud.nasa.gov/nsidc-cumulus-prod-protected/ATLAS/ATL03/006/2021/10/02/ATL03_20211002004127_01461306_006_01.h5"
#> [3,] "https://data.nsidc.earthdatacloud.nasa.gov/nsidc-cumulus-prod-protected/ATLAS/ATL03/006/2021/10/02/ATL03_20211002015115_01471302_006_01.h5"
#> [4,] "https://data.nsidc.earthdatacloud.nasa.gov/nsidc-cumulus-prod-protected/ATLAS/ATL03/006/2021/10/02/ATL03_20211002021545_01471306_006_01.h5"
#> [5,] "https://data.nsidc.earthdatacloud.nasa.gov/nsidc-cumulus-prod-protected/ATLAS/ATL03/006/2021/10/02/ATL03_20211002032533_01481302_006_01.h5"
#> [6,] "https://data.nsidc.earthdatacloud.nasa.gov/nsidc-cumulus-prod-protected/ATLAS/ATL03/006/2021/10/02/ATL03_20211002035002_01481306_006_01.h5"
```

Now we download the granules:


```r
# Download all granules
ATLAS_dataDownload(atl03_granules_local, "/the/directory/to/save")
```

And then we can open and work with them


```r
# Read the granules
atl03_h5 <- ATL03_read("/the/directory/to/save/name_of_granule.h5")

# List groups within atl03_h5
atl03_h5$ls()
```


```
#> [1] "gt1r"       "orbit_info"
```

## Working in the cloud

```r
atl03_granules_cloud <- ATLAS_dataFinder(
  short_name = "ATL03",
  lower_left_lon,
  lower_left_lat,
  upper_right_lon,
  upper_right_lat,
  version = "006",
  daterange = daterange,
  persist = TRUE,
  cloud_computing = TRUE
)

head(atl03_granules_cloud)
#> Collection: {'EntryTitle': 'ATLAS/ICESat-2 L2A Global Geolocated Photon Data V006'}
#> Spatial coverage: {'HorizontalSpatialDomain': {'Geometry': {'GPolygons': [{'Boundary': {'Points': [{'Longitude': 167.89473, 'Latitude': 59.54564}, {'Longitude': 167.67404, 'Latitude': 59.53425}, {'Longitude': 167.81095, 'Latitude': 58.84718}, {'Longitude': 168.29123, 'Latitude': 56.32957}, {'Longitude': 168.82548, 'Latitude': 53.24275}, {'Longitude': 169.43429, 'Latitude': 49.36326}, {'Longitude': 169.987, 'Latitude': 45.51084}, {'Longitude': 170.50035, 'Latitude': 41.65845}, {'Longitude': 170.98695, 'Latitude': 37.77973}, {'Longitude': 171.54559, 'Latitude': 33.0688}, {'Longitude': 172.15242, 'Latitude': 27.69247}, {'Longitude': 172.23362, 'Latitude': 26.9541}, {'Longitude': 172.35961, 'Latitude': 26.96509}, {'Longitude': 172.27919, 'Latitude': 27.70353}, {'Longitude': 171.67947, 'Latitude': 33.08007}, {'Longitude': 171.12888, 'Latitude': 37.7909}, {'Longitude': 170.65047, 'Latitude': 41.66957}, {'Longitude': 170.14702, 'Latitude': 45.52198}, {'Longitude': 169.60642, 'Latitude': 49.37448}, {'Longitude': 169.01275, 'Latitude': 53.25412}, {'Longitude': 168.49331, 'Latitude': 56.34106}, {'Longitude': 168.02745, 'Latitude': 58.85889}, {'Longitude': 167.89473, 'Latitude': 59.54564}]}}]}}}
#> Temporal coverage: {'RangeDateTime': {'BeginningDateTime': '2021-10-02T00:16:58.259Z', 'EndingDateTime': '2021-10-02T00:25:28.858Z'}}
#> Size(MB): 1732.5397911071777
#> Data: ['https://data.nsidc.earthdatacloud.nasa.gov/nsidc-cumulus-prod-protected/ATLAS/ATL03/006/2021/10/02/ATL03_20211002001658_01461302_006_01.h5']
```

In cloud computing you don't need to download data, instead you can 
read the data and start working with it.


```r
# Read the granule (the ATL03_read can only read one granule per read)
atl03_h5 <- ATL03_read(atl03_granules_cloud[1])
#> Opening 1 granules, approx size: 1.69 GB

# List groups within the h5 in cloud
atl03_h5$ls()
#>  [1] "METADATA"               "ancillary_data"         "atlas_impulse_response"
#>  [4] "ds_surf_type"           "ds_xyz"                 "gt1l"                  
#>  [7] "gt1r"                   "gt2l"                   "gt2r"                  
#> [10] "gt3l"                   "gt3r"                   "orbit_info"            
#> [13] "quality_assessment"
```


```r
# Which are strong beams
print(atl03_h5$strong_beams)
#> [1] "gt1l" "gt2l" "gt3l"
```


```r
# Orientation 0=forward, 1=backwards, 2=transition
print(atl03_h5[["orbit_info/sc_orient"]][])
#> [1] 1
```


```r
close(atl03_h5)
```












## Extract attributes


```r
# ATL03 seg attributes
atl03_seg_dt <- ATL03_seg_attributes_dt(atl03_h5, attributes = c("delta_time", "solar_elevation", "pitch", "h_ph", "ref_elev"))

head(atl03_seg_dt)
```



| delta_time| solar_elevation|      pitch| ref_elev| reference_photon_lon| reference_photon_lat|beam |
|----------:|---------------:|----------:|--------:|--------------------:|--------------------:|:----|
|  134086984|        33.53393| -0.0338702| 1.542873|            -106.5699|             41.53904|gt1r |
|  134086984|        33.53402| -0.0338686| 1.542873|            -106.5699|             41.53886|gt1r |
|  134086984|        33.53411| -0.0338672| 1.542873|            -106.5699|             41.53868|gt1r |
|  134086984|        33.53419| -0.0338657| 1.542873|            -106.5699|             41.53850|gt1r |
|  134086984|        33.53429| -0.0338641| 1.542873|            -106.5700|             41.53832|gt1r |
|  134086984|        33.53437| -0.0338627| 1.542873|            -106.5700|             41.53815|gt1r |



```r
# ATL08 seg attributes
atl08_seg_dt <- ATL08_seg_attributes_dt(atl08_h5, attributes = c("h_canopy", "h_te_mean", "terrain_slope", "canopy_openness"))

head(atl08_seg_dt)
```


| latitude| longitude|beam |  h_canopy| h_te_mean| terrain_slope| canopy_openness|
|--------:|---------:|:----|---------:|---------:|-------------:|---------------:|
| 41.53868| -106.5699|gt1r |  6.623291|  2448.891|    -0.0410579|       1.3961231|
| 41.53778| -106.5700|gt1r | 10.518555|  2447.029|     0.0256335|       2.0318608|
| 41.53689| -106.5701|gt1r |  6.695557|  2455.972|     0.0584960|       1.2851287|
| 41.53599| -106.5703|gt1r |  8.509766|  2462.762|     0.1706411|       2.2644644|
| 41.53509| -106.5704|gt1r |  4.614258|  2477.444|     0.0584462|       0.9932709|
| 41.53419| -106.5705|gt1r |  9.282227|  2485.570|     0.0930199|       2.1762936|


### Plot histograms:


```r
layout(t(1:2))

# ATL03 height histogram
hist(atl03_seg_dt$ref_elev, col = "#bd8421", xlab = "Elevation (m)", main = "h_te_mean")
hist(atl08_seg_dt$h_canopy, col = "green", xlab = "height (m)", main = "h_canopy")
```

<div align="center">

<div class="figure" style="text-align: center">
<img src="figure/unnamed-chunk-51-1.png" alt="Histograms for ATL03 elevation and ATL08 h_canopy"  />
<p class="caption">Histograms for ATL03 elevation and ATL08 h_canopy</p>
</div>

</div>

## Export to vector
The function `to_vect()` will return a `terra::vect` object.




```r
library(terra)

blueYellowRed <- function(n) grDevices::hcl.colors(n, "RdYlBu")

atl03_seg_vect <- to_vect(atl03_seg_dt)

# Plot with terra::plet
mapview::mapView(
  atl03_seg_vect,
  zcol = "ref_elev",
  breaks = 3,
  col.regions = blueYellowRed,
  map.types = c("Esri.WorldImagery")
)
```

<div align="center">


<img src="figure/atl03_seg_vect.png" width=500 />


</div>


```r
atl08_seg_vect <- to_vect(atl08_seg_dt)

greenYellowRed <- function(n) {
  grDevices::hcl.colors(n, "RdYlGn")
}

# Plot with terra::plet
map_vect <- mapview::mapView(
  atl08_seg_vect,
  zcol = "h_canopy",
  col.regions = greenYellowRed,
  map.types = c("Esri.WorldImagery")
)

map_vect
```

<div align="center">

<img src="figure/atl08_seg_vect.png" width=500 />



</div>

Save vector as geopackage file. The formats supported are as from GDAL terra package.


```r
terra::writeVector(atl03_seg_vect, "atl03_seg.gpkg")
terra::writeVector(atl08_seg_vect, "atl08_seg.gpkg")
```

## View ATL08 segments as raster

Single max_h_canopy:


```r
redYellowGreen <- function(n) grDevices::hcl.colors(n, "RdYlGn")
max_h_canopy <- ATL08_seg_attributes_dt_gridStat(atl08_seg_dt, func = max(h_canopy), res = 0.002)

mapview::mapView(
  max_h_canopy,
  map = map_vect,
  col.regions = redYellowGreen
)
```

<div align="center">

<img src="figure/atl08_max_h_canopy.png" width=500 />


</div>

### Multiple data:


```r
multiple_data <- ATL08_seg_attributes_dt_gridStat(atl08_seg_dt, func = list(
  max_h_canopy = max(h_canopy),
  min_h_canopy = min(h_canopy),
  mean_canopy_openness = mean(canopy_openness),
  mean_h_te_mean = mean(h_te_mean)
), res = 0.002)

map_vect_openness <- mapview::mapView(
  atl08_seg_vect,
  zcol = "canopy_openness",
  col.regions = redYellowGreen,
  map.types = c("Esri.WorldImagery")
)

blueYellowRed <- function(n) grDevices::hcl.colors(n, "RdYlBu", rev = TRUE)

map_vect_terrain <- mapview::mapView(
  atl08_seg_vect,
  zcol = "h_te_mean",
  col.regions = blueYellowRed,
  map.types = c("Esri.WorldImagery")
)


m1 <- mapview::mapView(multiple_data[[1]], map = map_vect, col.regions = redYellowGreen)
m2 <- mapview::mapView(multiple_data[[2]], map = map_vect, col.regions = redYellowGreen)
m3 <- mapview::mapView(multiple_data[[3]], map = map_vect_openness, col.regions = redYellowGreen)
m4 <- mapview::mapView(multiple_data[[4]], col.regions = blueYellowRed, map = map_vect_terrain)

leafsync::sync(m1, m2, m3, m4)
```

<div align="center" style="width:calc(100% - 20px);left:10px;position:absolute">







![multi](figure/output_multi.png)



</div>






## Extract attributes


```r
atl03_atl08_dt <- ATL03_ATL08_photons_attributes_dt_join(atl03_h5, atl08_h5, beam = "gt1r")

head(atl03_atl08_dt)
```



<div align="center" style="overflow-x: scroll;">


| ph_segment_id|    lon_ph|   lat_ph|     h_ph| quality_ph| solar_elevation| dist_ph_along| dist_ph_across| night_flag| classed_pc_indx|beam | classed_pc_flag|      ph_h| d_flag| delta_time|
|-------------:|---------:|--------:|--------:|----------:|---------------:|-------------:|--------------:|----------:|---------------:|:----|---------------:|---------:|------:|----------:|
|        771236| -106.5699| 41.53912| 2454.684|          0|        33.53393|      1.148141|       12579.52|          0|               6|gt1r |               2|  2.619385|      1|  134086984|
|        771236| -106.5699| 41.53911| 2453.936|          0|        33.53393|      1.854454|       12579.54|          0|              12|gt1r |               2|  2.018066|      1|  134086984|
|        771236| -106.5699| 41.53911| 2455.669|          0|        33.53393|      1.861945|       12579.50|          0|              13|gt1r |               2|  3.890869|      1|  134086984|
|        771236| -106.5699| 41.53911| 2452.849|          0|        33.53393|      2.560772|       12579.57|          0|              20|gt1r |               2|  1.206787|      1|  134086984|
|        771236| -106.5699| 41.53909| 2459.739|          0|        33.53393|      4.717131|       12579.39|          0|              45|gt1r |               3|  8.225098|      1|  134086984|
|        771236| -106.5699| 41.53909| 2449.157|          0|        33.53393|      4.676751|       12579.68|          0|              46|gt1r |               0| -2.237305|      1|  134086984|

</div>

## Plotting the result:

```r
oldpar <- par(no.readonly = TRUE)
par(oma = c(0, 0, 0, 0))
par(mar = c(2, 3, 1, 1))
layout(matrix(c(1, 2), ncol = 1))
plot(
  atl03_atl08_dt,
  y = "h_ph",
  colors = c("gray", "#bd8421", "forestgreen", "green"),
  xlim = range(atl03_atl08_dt$dist_ph_along),
  ylim = range(atl03_atl08_dt$h_ph),
  beam = "gt1r",
  cex = 0.5,
  pch = 16
)

par(mar = c(3, 3, 1, 1))

plot(
  atl03_atl08_dt,
  y = "ph_h",
  colors = c("gray", "#bd8421", "forestgreen", "green"),
  xlim = range(atl03_atl08_dt$dist_ph_along),
  ylim = range(atl03_atl08_dt$ph_h),
  beam = "gt1r",
  cex = 0.5,
  pch = 16
)

par(
  oldpar
)
```

<div align="center">

<div class="figure" style="text-align: center">
<img src="figure/unnamed-chunk-96-1.png" alt="Classified ATL03 photons using ATL08 labels"  />
<p class="caption">Classified ATL03 photons using ATL08 labels</p>
</div>

</div>

## Calculating raster statistics


```r
h_canopy <- ATL03_ATL08_photons_attributes_dt_gridStat(
  atl03_atl08_dt,
  func = list(h_canopy = quantile(ph_h, 0.98), count = .N),
  res = 0.0001
)

plot(h_canopy,
  col = viridis::inferno(100),
  xlab = "Langitude (degree)",
  ylab = "Latitude (degree)"
)
```

<div align="center">

<div class="figure" style="text-align: center">
<img src="figure/unnamed-chunk-98-1.png" alt="Rasterized ATL03_ATL08 data for canopy height (h_canopy) and number of photons (n)"  />
<p class="caption">Rasterized ATL03_ATL08 data for canopy height (h_canopy) and number of photons (n)</p>
</div>

</div>









In this example we will model the `h_canopy` of the ICESat-2 using only the Harmonized Landsat Sentinel-2 dataset (hls).


## Extract ATL08 segment attributes h_canopy attribute


```r
atl08_seg_dt <- ATL08_seg_attributes_dt(atl08_h5, attribute = "h_canopy")
head(atl08_seg_dt)
```




```r
head(atl08_seg_dt)
```



| latitude| longitude|beam |  h_canopy|
|--------:|---------:|:----|---------:|
| 41.53868| -106.5699|gt1r |  6.623291|
| 41.53778| -106.5700|gt1r | 10.518555|
| 41.53689| -106.5701|gt1r |  6.695557|
| 41.53599| -106.5703|gt1r |  8.509766|
| 41.53509| -106.5704|gt1r |  4.614258|
| 41.53419| -106.5705|gt1r |  9.282227|



### Visualizing the 'h_canopy' for the ATL08 dataset.


```r
library(terra)

atl08_seg_vect <- to_vect(atl08_seg_dt)
terra::plet(atl08_seg_vect, "h_canopy", col = grDevices::hcl.colors(9, "RdYlGn"), tiles = c("Esri.WorldImagery"))
```

<div align="center" style="display:flex;justify-content:center">



<img src="figure/atl08_seg_vect_gee_modelling.png" width=500 />


</div>

### Querying the GEEs datasets for Harmonized Landsat Sentinel-2


```r
hls_search <- search_datasets("Harmonized", "Landsat")
hls_search
#>                      id
#>                  <char>
#> 1: NASA_HLS_HLSL30_v002
#>                                                                                                    title
#>                                                                                                   <char>
#> 1: HLSL30: HLS-2 Landsat Operational Land Imager Surface Reflectance and TOA Brightness Daily Global 30m
#>                                                                                                                                                                                                                                                                                                                          description
#>                                                                                                                                                                                                                                                                                                                               <char>
#> 1: The Harmonized Landsat Sentinel-2 (HLS) project provides consistent surface reflectance (SR) and top of atmosphere (TOA) brightness data from a virtual constellation of satellite sensors. The Operational Land Imager (OLI) is housed aboard the joint NASA/USGS Landsat 8 and Landsat 9 satellites, while the Multi-Spectral …
```


```r
hls_id <- get_catalog_id(hls_search$id)
hls_id
#> [1] "NASA/HLS/HLSL30/v002"
```

### Open the Google Earth Engine HLS catalog and get band names


```r
hls_collection <- ee$ImageCollection(hls_id)
names(hls_collection)
#>  [1] "B1"    "B2"    "B3"    "B4"    "B5"    "B6"    "B7"    "B9"    "B10"  
#> [10] "B11"   "Fmask" "SZA"   "SAA"   "VZA"   "VAA"
```

### Define area of interest (aoi) clip boundaries and time and cloud mask for filtering.


```r
bbox <- terra::ext(atl08_seg_vect)

aoi <- ee$Geometry$BBox(
  west = bbox$xmin,
  south = bbox$ymin,
  east = bbox$xmax,
  north = bbox$ymax
)

hls <- hls_collection$
  filterDate("2022-04-01", "2022-04-30")$
  filterBounds(aoi)$
  map(function(x) x$updateMask(!(x[["Fmask"]] & 14)))$
  median()


hls_unmasked <- hls_collection$
  filterDate("2022-04-01", "2022-04-30")$
  filterBounds(aoi)$
  median()
```

### Calculate EVI:


```r
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
#> ee.image.Image
#> 
#> Bands
#> [1] "blue"  "green" "red"   "nir"   "swir1" "swir2" "evi"
```

## Visualize the resulting image


```r
library(leaflet)

forest_height_palette <- c("#ffffff", "#8b4513", "#99cc99", "#006600", "#004d00")
palette_colors <- colorNumeric(forest_height_palette, range(atl08_seg_dt$h_canopy))(atl08_seg_dt[order(h_canopy), h_canopy])

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
```

<div align="center" style="display:flex;justify-content:center">

<img src="figure/mask_unmasked_hls.png" width=500 />


</div>

## Extracting GEE data for segments

For each segment extract the hls data:


```r
extracted_dt <- seg_gee_ancillary_dt_extract(hls, atl08_seg_vect)
#> Processing 1-9 of 9

head(extracted_dt)
```



| idx|beam |  h_canopy|     red|   green|    blue|     nir|   swir1|   swir2|        evi|
|---:|:----|---------:|-------:|-------:|-------:|-------:|-------:|-------:|----------:|
|   1|gt1r |  6.623291| 0.69625| 0.67005| 0.64725| 0.67485| 0.05005| 0.03855| -0.0536086|
|   2|gt1r | 10.518555| 0.43980| 0.42210| 0.39670| 0.51935| 0.05140| 0.03460|  0.1681249|
|   3|gt1r |  6.695557| 0.39235| 0.39130| 0.36095| 0.47015| 0.04505| 0.02975|  0.1741076|
|   4|gt1r |  8.509766| 0.34855| 0.34120| 0.31205| 0.40995| 0.05280| 0.03950|  0.1322278|
|   5|gt1r |  4.614258| 0.47720| 0.46660| 0.45450| 0.49330| 0.04385| 0.03400|  0.0424690|
|   6|gt1r |  9.282227| 0.24345| 0.24430| 0.20735| 0.36385| 0.06070| 0.04180|  0.2371152|




## Fit the randomForest model


```r
bandNames <- names(hls)
x <- extracted_dt[, .SD, .SDcols = bandNames]
y <- extracted_dt[["h_canopy"]]

rf_model <- model_fit(x, y, ntree = 500, mtry = 1)
print(rf_model)
#> 
#> Call:
#>  randomForest(x = x, y = y, ntree = 500, mtry = 1) 
#>                Type of random forest: regression
#>                      Number of trees: 500
#> No. of variables tried at each split: 1
#> 
#>           Mean of squared residuals: 4.206443
#>                     % Var explained: -57.46
```


```r
library(randomForest)

rf_importance <- importance(rf_model)
barplot(rf_importance[, "IncNodePurity"], ylim = c(0, 3), main = "Variable importance (Increase Node Purity)")
```

<div align="center">

<div class="figure" style="text-align: center">
<img src="figure/unnamed-chunk-129-1.png" alt="Random forests variable importance (increase node impurity)." width="500" />
<p class="caption">Random forests variable importance (increase node impurity).</p>
</div>

</div>

## Apply the model to Google Earth Engine WorldImagery


```r
gee_model <- build_ee_forest(rf_model)
result <- hls$classify(gee_model)
min_hcanopy <- min(atl08_seg_dt$h_canopy)
max_hcanopy <- max(atl08_seg_dt$h_canopy)
atl08_seg_vect$h_canopy <- round(atl08_seg_vect$h_canopy, 3) # Round off to 3 decimal places

modelled_map <- terra::plet(
  atl08_seg_vect,
  "h_canopy",
  palette_colors,
  tiles = ""
) |>
  addEEImage(
    hls,
    bands = c("red", "green", "blue"),
    group = "hls",
    min = 0,
    max = 0.6
  ) |>
  addEEImage(
    result,
    bands = "classification",
    group = "classification",
    min = min_hcanopy,
    max = max_hcanopy,
    palette = forest_height_palette
  ) |>
  leaflet::addLegend(
    pal = colorNumeric(forest_height_palette, seq(min_hcanopy, max_hcanopy)),
    values = seq(min_hcanopy, max_hcanopy, length = 3),
    opacity = 1,
    title = "h_canopy",
    position = "bottomleft",
  ) |>
  setView(lng = centroid[1], lat = centroid[2], zoom = 13) |>
  addLayersControl(
    overlayGroups = c("classification"),
    options = layersControlOptions(collapsed = FALSE)
  )

modelled_map
```

<div align="center" style="display:flex;justify-content:center">

<img src="figure/upscalled_gee_map.png" width=500 />


</div>



# Acknowledgements
We gratefully acknowledge funding from NASA’s ICESat-2 (ICESat-2, grant 22-ICESat2_22-0006), Carbon Monitoring System (CMS, grant 22-CMS22-0015) and Commercial Smallsat Data Scientific Analysis(CSDSA, grant 22-CSDSA22_2-0080). 

# Reporting Issues 
Please report any issue regarding the ICESat2VegR package to Dr. Silva (c.silva@ufl.edu)

# Citing ICESat2VegR
Silva,C.A; Hamamura,C.ICESat2VegR: An R Package for NASA's Ice, Cloud, and Elevation Satellite (ICESat-2) Data Processing and Visualization for Terrestrial Applications.version 0.0.1, accessed on November. 22 2023, available at: <https://CRAN.R-project.org/package=ICESat2VegR>

# Disclaimer
**ICESat2VegR package comes with no guarantee, expressed or implied, and the authors hold no responsibility for its use or reliability of its outputs.**

