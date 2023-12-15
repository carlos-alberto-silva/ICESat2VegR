![](https://github.com/carlos-alberto-silva/rICESat2Veg/blob/master/readme/cover.png)<br/>
[![R-CMD-check](https://github.com/carlos-alberto-silva/rICESat2Veg/actions/workflows/r.yml/badge.svg?branch=master)](https://github.com/carlos-alberto-silva/rICESat2Veg/actions/workflows/r.yml)
[![CRAN](https://www.r-pkg.org/badges/version/rICESat2Veg)](https://cran.r-project.org/package=rICESat2Veg)
![Github](https://img.shields.io/badge/Github-0.1.12-green.svg)
![licence](https://img.shields.io/badge/Licence-GPL--3-blue.svg) 
![Downloads](https://cranlogs.r-pkg.org/badges/grand-total/rICESat2Veg)
[![Build Status](https://travis-ci.com/carlos-alberto-silva/rICESat2Veg.svg?token=Jqizwyc6gBxNafNccTdU&branch=master)](https://travis-ci.com/carlos-alberto-silva/rICESat2Veg)

**rICESat2Veg: An R Package for NASA's Ice, Cloud, and Elevation Satellite (ICESat-2) Data Processing and Visualization for Land and Vegetation Applications.**

Authors: Carlos Alberto Silva and Caio Hamamura  

The rICESat2Veg package provides functions for downloading, processing, visualizing and exporting ICESat-2 ATL03 and ATL08 data in R environment.

# Getting Started

## Installation
```r
#The CRAN version:
install.packages("rICESat2Veg")

# The development version:
#install.packages("remotes")
remotes::install_github("https://github.com/carlos-alberto-silva/rICESat2Veg", dependencies = TRUE)

# loading rGEDI package
library(rICESat2Veg)

```    
## Downloading ICESat-2 data
```r
#######
# Herein, we are using only a ICESat-2 sample dataset for this tutorial.
#######
# downloading zip file
download.file("https://github.com/carlos-alberto-silva/rICESat2Veg/releases/download/datasets/examples.zip",destfile=file.path(outdir, "examples.zip"))

# unzip file 
unzip(file.path(outdir,"examples.zip"))

```

## Reading ATL03 and ATL08 h5 files
```r
atl03_h5<-ATL03_read(atl03_path=atl03_path)
atl08_h5<-ATL08_read(atl08_path=atl08_path)
```

## Extracting ATL03 and ATL08 photon attributes
```r
atl03_photons_dt<-ATL03_photons_attributes_dt(atl03_h5=atl03_h5, beam="gt1r")
head(atl03_photons_dt)

##     lon_ph   lat_ph     h_ph     quality_ph solar_elevation dist_ph_along
## 1: -103.7555 59.49127 422.8559          0        23.49145     0.5788960
## 2: -103.7555 59.49127 413.9451          0        23.49145     0.5501707
## 3: -103.7555 59.49127 301.0684          0        23.49145     0.1823875
## 4: -103.7555 59.49126 434.1847          0        23.49145     1.3294470
## 5: -103.7555 59.49126 398.2455          0        23.49145     1.2123015
## 6: -103.7555 59.49126 369.2183          0        23.49145     1.1181128

atl08_photons_dt<-ATL08_photons_attributes_dt(atl08_h5=atl08_h5, beam="gt1r")
head(atl08_photons_dt)

##    ph_segment_id beam classed_pc_indx classed_pc_flag  ph_h     d_flag   delta_time
## 1: 671089        gt1r               8               1 -0.1245117      1  134086702
## 2: 671089        gt1r               9               1 -0.2919922      1  134086702
## 3: 671089        gt1r              11               0 -1.0151367      1  134086702
## 4: 671089        gt1r              15               0 -0.8881836      1  134086702
## 5: 671089        gt1r              23               0 -0.9491882      1  134086702
## 6: 671089        gt1r              33               0 -0.6138306      1  134086702

```

## Extracting ATL08-derived terrain and canopy attributes by segments
```r
# Extracting
atl08_seg_att_dt <- ATL08_seg_attributes_dt(atl08_h5 = atl08_h5)
head(atl08_seg_att_dt)

##   latitude longitude beam h_canopy canopy_openness h_te_mean terrain_slope
## 1: 59.49081 -103.7556 gt1r 4.820282       1.0603539  343.2197    0.06025224
## 2: 59.48993 -103.7558 gt1r 2.681976       0.7612819  345.8327   -0.01439811
## 3: 59.48903 -103.7560 gt1r 3.851379       0.9414949  346.6401    0.04266450
## 4: 59.48814 -103.7561 gt1r 4.076965       0.9898106  350.9311    0.04225401
## 5: 59.48725 -103.7563 gt1r 6.157623       1.1411711  350.4123   -0.11771229
## 6: 59.48635 -103.7565 gt1r 6.686523       1.4989151  340.1538   -0.05887743


summary(atl08_seg_att_dt)

# plotting histograms
atl08_seg_att_dt@dt$h_canopy[atl08_seg_att_dt@dt$h_canopy>80]<-NA # set NA to values > 80 m
atl08_seg_att_dt@dt$h_te_mean[atl08_seg_att_dt@dt$h_te_mean>5000]<-NA # set NA to values > 5000 m
par(mfrow=c(1,2))
hist(atl08_seg_att_dt$h_canopy, col="green", xlab="height (m)", main="h_canopy")
hist(atl08_seg_att_dt$h_te_mean, col="#bd8421", xlab="Elevation (m)", main="h_te_mean")
```
![](https://github.com/carlos-alberto-silva/rICESat2Veg/blob/master/readme/hist_ATL08.png)

```r
# Plotting ATL08 attribute on a map by segment coordinates
library(leaflet) # loading leaflet package
library(leafsync) # loading leafsync package

# Set breaks values for h_canopy and h_te_mean
options(scipen=999)
h_canopy.breaks=seq(0,60,10)
h_te_mean.breaks=seq(3000,4000,200)

h_canopy_bins<-h_canopy.breaks[cut(atl08_seg_att_dt$h_canopy,
                   breaks=h_canopy.breaks,labels=FALSE)]

h_te_mean_bins<-h_te_mean.breaks[cut(atl08_seg_att_dt$h_te_mean,
                          breaks=h_te_mean.breaks,labels=FALSE)]

# set color palette
h_canopy.pal <- colorBin(palette = "BrBG", domain = h_canopy_bins)
h_te_mean.pal <- colorBin(palette = "BrBG", domain = h_te_mean_bins)

# plotting h_canopy
h_canopy.map<-leaflet() %>%
     setView(lng = -107.11, lat = 37.51314,zoom = 13.49) %>%
     addCircleMarkers(atl08_seg_att_dt$longitude,
                      atl08_seg_att_dt$latitude,
       radius = 0.5,
       opacity = 1,
       color = h_canopy.pal(h_canopy_bins)
     ) %>%
     addScaleBar(options = list(imperial = FALSE)) %>%
     addProviderTiles(providers$Esri.WorldImagery) %>%
     addLegend(pal = h_canopy.pal,values = h_canopy_bins,title ="h_canopy (m)")

# plotting h_te_mean
h_te_mean.map<-leaflet() %>%
  setView(lng = -107.11, lat = 37.51314,zoom = 13.49) %>%
  addCircleMarkers(atl08_seg_att_dt$longitude,
                   atl08_seg_att_dt$latitude,
                   radius = 0.5,
                   opacity = 1,
                   color = h_te_mean.pal(h_te_mean_bins)
  ) %>%
  addScaleBar(options = list(imperial = FALSE)) %>%
  addProviderTiles(providers$Esri.WorldImagery) %>%
  addLegend(pal = h_te_mean.pal,values = h_te_mean_bins,title ="h_te_mean (m)")
sync(h_canopy.map, h_te_mean.map)
```
![](https://github.com/carlos-alberto-silva/rICESat2Veg/blob/master/readme/ATL08_segments.png)

## Clipping ATL08 Terrain and Canopy Attributes

```r
# Clipping by  bounding box
# Define the bounding box
xmin <- -107.7
xmax <- -106.5
ymin <- 32.75
ymax <- 42.75

# Clip
atl08_seg_att_dt_clip <- ATL08_seg_attributes_dt_clipBox(atl08_seg_att_dt, xmin, xmax, ymin, ymax)
head(atl08_seg_att_dt_clip) # print the first six observations

##   latitude longitude beam  h_canopy canopy_openness h_te_mean  terrain_slope
## 1: 42.08012 -106.5000 gt1r  4.976807        1.069874  2075.281 -0.11553377658
## 2: 42.07922 -106.5001 gt1r  7.460938        1.273932  2063.163 -0.06796860695
## 3: 42.07832 -106.5003 gt1r 10.185303        2.310288  2056.079 -0.08571143448
## 4: 42.07742 -106.5004 gt1r  8.583984        2.444624  2050.611 -0.02976966091
## 5: 42.07652 -106.5005 gt1r  6.843018        2.086431  2048.529 -0.00003621415
## 6: 42.07563 -106.5006 gt1r  6.016113        1.481141  2049.141  0.00491552008

# Clipping by geometry
# Specify the path to shapefile
poly_filepath <- system.file("extdata", "polygon.shp", package = "rICESat2Veg")

# Read shapefile
poly <- terra::vect(poly_filepath)


# Clip
atl08_seg_att_dt_clipg <- ATL08_seg_attributes_dt_clipGeometry(atl08_seg_att_dt, poly, split_by = "FID")
head(atl08_seg_att_dt_clipg) # print the first six observations

##   latitude longitude beam h_canopy canopy_openness h_te_mean terrain_slope nid  poly_id
## 1: 52.49692 -105.0109 gt1r 4.191498        1.060516  503.8196  0.0012361892   1       1
## 2: 52.47792 -105.0140 gt1r 7.399384        1.957333  513.7219  0.1044191048   2       1
## 3: 52.47721 -105.0141 gt1r 3.816345        1.160930  514.6264  0.0050510662   3       1
## 4: 52.46538 -105.0160 gt1r 6.513763        1.305076  504.8456 -0.0748032928   4       1
## 5: 52.46466 -105.0161 gt1r 6.364258        1.470707  503.5898  0.0159635153   5       1
## 6: 52.46377 -105.0163 gt1r 4.756866        1.288144  503.8241 -0.0004335443   6       1

## View ATL08 clipped data by bbox
m1<-leaflet() %>%
  setView(lng = -107.11, lat = 37.51314,zoom = 05) %>%
  addCircleMarkers(atl08_seg_att_dt$longitude,
                   atl08_seg_att_dt$latitude,
                   radius = 1,
                   opacity = 1,
                   color = "red")  %>%
  addCircleMarkers(atl08_seg_att_dt_clip$longitude,
                   atl08_seg_att_dt_clip$latitude,
                   radius = 1,
                   opacity = 1,
                   color = "green")  %>%
  addScaleBar(options = list(imperial = FALSE)) %>%
  addProviderTiles(providers$Esri.WorldImagery)  %>%
  addLegend(colors = c("red","green"), labels= c("All samples","Clip bbox"),title ="ATL08")

## View ATL08 clipped data by geometry
# color palette
pal <- colorFactor(
  palette = c( 'green','orange'),
  domain = atl08_seg_att_dt_clipg$poly_id
)

m2<-leaflet() %>%
  setView(lng = -107.11, lat = 37.51314,zoom = 03) %>%
  addCircleMarkers(atl08_seg_att_dt$longitude,
                   atl08_seg_att_dt$latitude,
                   radius = 1,
                   opacity = 1,
                   color = "red")  %>%
  addCircleMarkers(atl08_seg_att_dt_clipg$longitude,
                   atl08_seg_att_dt_clipg$latitude,
                   radius = 1,
                   opacity = 1,
                   color = pal(atl08_seg_att_dt_clipg$poly_id))  %>%
  addScaleBar(options = list(imperial = FALSE)) %>%
  addPolygons(data=poly,weight=1,col = 'white',
              opacity = 1, fillOpacity = 0) %>%
  addProviderTiles(providers$Esri.WorldImagery) %>%
  addLegend(pal = pal, values = atl08_seg_att_dt_clipg$poly_id,title ="Poly IDs" )
sync(m1, m2)
```
![](https://github.com/carlos-alberto-silva/rICESat2Veg/blob/master/readme/Fig_clipping_ATL08.png)

## Gridding AT08 Terrain and Canopy Attributes
```r
# Computing the max h_canopy at 0.05 degree grid cell
max_h_canopy <-ATL08_seg_attributes_dt_gridStat(atl08_seg_att_dt, func=max(h_canopy), res=0.05)
plot(max_h_canopy, xlim=c(-107.2,-106.8),ylim=c(38,39), col=viridis::inferno(8),
     main="Max 'h_canopy'",
     xlab="Langitude (degree)",
     ylab="Latitude (degree)")
```
<img align="right" src="https://github.com/carlos-alberto-silva/rICESat2Veg/blob/master/readme/grid_h_canopy.png"  width="300">

```r
# Define your own function
mySetOfMetrics <- function(x) {
  metrics <- list(
    min = min(x), # Min of x
    max = max(x), # Max of x
    mean = mean(x), # Mean of x
    sd = sd(x) # Sd of x
  )
  return(metrics)
}

# Computing a series of h_canopy statistics at 0.05 degree grid cell from customized function
atl08_seg_att_dt@dt$h_canopy[atl08_seg_att_dt@dt$h_canopy>80]<-NA # set values > 80 m to NA m
h_canopy_metrics <-ATL08_seg_attributes_dt_gridStat(atl08_seg_att_dt, func=mySetOfMetrics(h_canopy),res=0.05)
plot(h_canopy_metrics,
     xlim=c(-107.2,-106.8),
     ylim=c(38,39),
     col=viridis::viridis(8),
     xlab="Langitude (degree)",
     ylab="Latitude (degree)")
```
<img align="right" src="https://github.com/carlos-alberto-silva/rICESat2Veg/blob/master/readme/grid_h_canopy2.png"  width="300">


## Join ATL03 and ATL08 photon attributes tables
```r
atl03_atl08_dt<-ATL03_ATL08_photons_attributes_dt_join(atl03_h5,atl08_h5, beam = "gt1l")
head(atl03_atl08_dt)

## lon_ph   lat_ph     h_ph quality_ph solar_elevation dist_ph_along ph_segment_id classed_pc_indx
## 1: -103.7541 59.49119 375.0020          0        23.49148     0.7613657            NA              NA
## 2: -103.7541 59.49119 389.0057          0        23.49148     0.7323560            NA              NA
## 3: -103.7540 59.49118 339.6367          0        23.49148     0.8367972        671089               3
## 4: -103.7540 59.49118 333.1566          0        23.49148     0.8507242            NA              NA
## 5: -103.7540 59.49118 339.5233          0        23.49148     0.8367975        671089               5
## 6: -103.7540 59.49118 339.7158          0        23.49148     0.8367971        671089               6
## classed_pc_flag d_flag delta_time        ph_h beam night_flag
## 1:              NA     NA         NA          NA gt1l          1
## 2:              NA     NA         NA          NA gt1l          1
## 3:               1      1  134086701 -0.04504395 gt1l          1
## 4:              NA     NA         NA          NA gt1l          1
## 5:               1      1  134086701 -0.16015625 gt1l          1
## 6:               1      1  134086701  0.03042603 gt1l          1

## Plotting photon cloud
# plot by "ph_h"
plot(atl03_atl08_dt, y="ph_h",colors, beam = "gt1l",
     ylim=c(0,7), xlim=c(2400,3800), pch=16)
```
![](https://github.com/carlos-alberto-silva/rICESat2Veg/blob/master/readme/cloud_height.png)

```r
# plot by "h_ph"
plot(atl03_atl08_dt, y="h_ph",colors, beam = "gt1l",
     ylim=c(300,355),xlim=c(2400,3800), pch=16)
```
![](https://github.com/carlos-alberto-silva/rICESat2Veg/blob/master/readme/cloud_height2.png)

## Gridding ATL03 and AT08 Photons Attributes
```r
# Computing the mean of h_ph attribute at 0.05 degree grid cell
mean_h_ph <- ATL03_ATL08_photons_attributes_dt_gridStat(atl03_atl08_dt, func = mean(h_ph), res = 0.05)
plot(mean_h_ph,
     xlim=c(-107.2,-106.8),
     ylim=c(38,39),
     main="Mean h_ph",
     col=viridis::inferno(8),
     xlab="Langitude (degree)",
     ylab="Latitude (degree)")
```
<img align="right" src="https://github.com/carlos-alberto-silva/rICESat2Veg/blob/master/readme/h_ph_map.png"  width="300">

```r
# Define your own function
mySetOfMetrics <- function(x) {
  metrics <- list(
    min = min(x), # Min of x
    max = max(x), # Max of x
    mean = mean(x), # Mean of x
    sd = sd(x) # Sd of x
  )
  return(metrics)
}

# Computing a series of h_canopy statistics at 0.05 degree grid cell from customized function
h_ph_metrics <-ATL03_ATL08_photons_attributes_dt_gridStat(atl03_atl08_dt, func=mySetOfMetrics(h_ph),res=0.05)
plot(h_ph_metrics,
     xlim=c(-107.2,-106.8),
     ylim=c(38,39),
     col=viridis::viridis(8),
     xlab="Langitude (degree)",
     ylab="Latitude (degree)")
```
<img align="right" src="https://github.com/carlos-alberto-silva/rICESat2Veg/blob/master/readme/h_ph_map_all.png"  width="300">

## Clipping ATL03 and joined ATL03 and ATL08 photons attributes

```r
# Clipping by  bounding box

# ATL03 photons
atl03_photons_dt_clip <- ATL03_photons_attributes_dt_clipBox(atl03_photons_dt, xmin, xmax, ymin, ymax)
head(atl03_photons_dt_clip) # print the first six observations
##   lon_ph   lat_ph     h_ph quality_ph solar_elevation dist_ph_along
## 1: -106.5 42.07277 2277.565          0        33.26923       1946643
## 2: -106.5 42.07275 2270.703          0        33.26932       1946645
## 3: -106.5 42.07274 2223.677          0        33.26932       1946647
## 4: -106.5 42.07273 2262.953          0        33.26932       1946648
## 5: -106.5 42.07271 2280.661          0        33.26932       1946649
## 6: -106.5 42.07270 2271.058          0        33.26932       1946652

# AT03 and ATL08 photons
atl03_atl08_photons_dt_clip <- ATL03_ATL08_photons_attributes_dt_clipBox(atl03_atl08_dt, xmin, xmax, ymin, ymax)
head(atl03_atl08_dt_clip) # print the first six observations
##      lon_ph   lat_ph     h_ph quality_ph solar_elevation dist_ph_along ph_segment_id
## 1: -106.5 42.07237 2082.053          0         33.2695       1946688        768267
## 2: -106.5 42.07236 2085.700          0         33.2695       1946688        768267
## 3: -106.5 42.07235 2091.462          0         33.2695       1946689        768267
## 4: -106.5 42.07235 2088.416          0         33.2695       1946689        768267
## 5: -106.5 42.07233 2081.228          0         33.2695       1946691        768267
## 6: -106.5 42.07233 2084.664          0         33.2695       1946691        768267
   classed_pc_indx classed_pc_flag d_flag delta_time       ph_h beam night_flag
## 1:              24               1      1  134086975 -0.2282715 gt1l          1
## 2:              28               2      1  134086975  3.2006836 gt1l          1
## 3:              33               3      1  134086975  8.7456055 gt1l          1
## 4:              34               3      1  134086975  5.4833984 gt1l          1
## 5:              44               0      1  134086975 -1.9194336 gt1l          1
## 6:              46               2      1  134086975  1.3044434 gt1l          1


# Clipping by geometry

# ATL03 photons
atl03_photons_dt_clipg <- ATL03_photons_attributes_dt_clipGeometry(atl03_photons_dt, polygon, split_by = "FID")
head(atl03_photons_dt_clipg) # print the first six observations

# ATL03 and ATL08 photons
atl03_atl08_photons_dt_clipg <- ATL03_ATL08_photons_attributes_dt_cliGeometry(atl03_photons_dt, polygon, split_by = "FID")
head(atl03_photons_dt_clipg) # print the first six observations
```

## Computing terrain and canopy metrics by a within segment length

```r
# Computing the max canopy height at 30 m segments
max_ph_h_seg_30m <-ATL03_ATL08_compute_seg_attributes_dt_segStat(atl03_atl08_dt, func=max(ph_h),
                                seg_length = 30,
                                ph_class=c(2,3),
                                beam=c("gt1l", "gt1r", "gt2l", "gt2r", "gt3l", "gt3r"),
                                quality_ph=0,
                                night_flag=1)
head(max_ph_h_seg_30m)
##      segment max_ph_h latitude longitude dist_along
## 1:       1 2.025024 59.49107 -103.7541         15
## 2:       2 2.578125 59.49079 -103.7541         45
## 3:       3 3.153015 59.49051 -103.7542         75
## 4:       4 2.994751 59.49026 -103.7542        105
## 5:       5 2.130676 59.49000 -103.7543        135
## 6:       6 2.302063 59.48966 -103.7544        165

# Computing the max canopy height at 100 m segments
max_ph_h_seg_100m <-ATL03_ATL08_compute_seg_attributes_dt_segStat(atl03_atl08_dt, func=max(ph_h),
                                                                 seg_length = 100,
                                                                 ph_class=c(2,3),
                                                                 beam=c("gt1l", "gt1r", "gt2l", "gt2r", "gt3l", "gt3r"),
                                                                 quality_ph=0,
                                                                 night_flag=1)

head(max_ph_h_seg)
##   segment max_ph_h latitude longitude dist_along
## 1:   1 2.025024 59.49107 -103.7541         15
## 2:   2 2.578125 59.49079 -103.7541         45
## 3:   3 3.153015 59.49051 -103.7542         75
## 4:   4 2.994751 59.49026 -103.7542        105
## 5:   5 2.130676 59.49000 -103.7543        135
## 6:   6 2.302063 59.48966 -103.7544        165


## Plot photon cloud and computed metrics
colors <- c("gray", "#bd8421", "forestgreen", "green")
colorMap <- colors[atl03_atl08_dt$classed_pc_flag + 1]
plot(atl03_atl08_dt$dist_ph_along[20000:30000],
     atl03_atl08_dt$ph_h[20000:30000],
     ylim=c(-2,8),
     col=colorMap[20000:30000],
     pch=16,
     xlab="Distance along-track (m)", ylab="Height (m)")
grid()
points(max_ph_h_seg_30m$dist_along,max_ph_h_seg_30m$max_ph_h, col="red", pch=16)
points(max_ph_h_seg_100m$dist_along,max_ph_h_seg_100m$max_ph_h, col="blue", pch=16)
legend("topleft",legend=c("ATL03 unclassified","ATL03 ground","ATL03 Canopy","ATL03 Top canopy","Max ph_h at 30m","Max ph_h at 100m"), col=c(colors,"red","blue"), pch=16, bty="n")

# Define your own function
 mySetOfMetrics <- function(x) {
   metrics <- list(
     min = min(x), # Min of x
     max = max(x), # Max of x
     mean = mean(x), # Mean of x
     sd = sd(x) # Sd of x
     h_canopy = quantile(x,0.98)
   )
   return(metrics)
}

# Computing a series of ph_h statistics from customized function
 ph_h_metrics <-ATL03_ATL08_compute_seg_attributes_dt_segStat(atl03_atl08_dt, func=mySetOfMetrics(ph_h),
                                seg_length = 30, # 30 m
                                ph_class=c(2,3),
                                beam=c("gt1l", "gt1r", "gt2l", "gt2r", "gt3l", "gt3r"),
                                quality_ph=0,
                                night_flag=1)

head(ph_h_metrics)
##  segment  min_ph_h max_ph_h mean_ph_h   sd_ph_h latitude longitude dist_along
## 1:  1 0.5019531 2.025024 0.9590389 0.5434316 59.49107 -103.7541         15
## 2:  2 0.5261841 2.578125 0.9932157 0.5950008 59.49079 -103.7541         45
## 3:  3 0.5100403 3.153015 1.1455467 0.7921570 59.49051 -103.7542         75
## 4:  4 0.5131226 2.994751 1.3727740 0.8207025 59.49026 -103.7542        105
## 5:  5 0.5229187 2.130676 1.1009442 0.4158372 59.49000 -103.7543        135
## 6:  6 0.5034790 2.302063 1.2780537 0.4787862 59.48966 -103.7544        165

```


## Statistics of ATL08 Terrain and Canopy Attributes by Geometry
```r
# Computing the max h_canopy by polygon id
max_h_canopy <-ATL08_seg_attributes_dt_polyStat(atl08_seg_att_dt_clipg, func=max(h_canopy),poly_id="poly_id")
head(max_h_canopy)
##    poly_id  max
## 1:  1      74.25031
## 2:  0      71.87158

# Computing a series of canopy statistics from customized function
h_canopy_metrics <-ATL08_seg_attributes_dt_polyStat(atl08_seg_att_dt_clipg, func=mySetOfMetrics(h_canopy),poly_id="poly_id")
head(h_canopy_metrics)
##   poly_id       min      max      mean       sd
## 1:       1 0.9290771 74.25031  7.860409 4.955381
## 2:       0 1.4843750 71.87158 10.627854 5.201939
```

# Acknowledgements
We gratefully acknowledge funding from NASA’s ICESat-2 (ICESat-2, grant 22-ICESat2_22-0006), Carbon Monitoring System (CMS, grant 22-CMS22-0015) and Commercial Smallsat Data Scientific Analysis(CSDSA, grant 22-CSDSA22_2-0080). 

# Reporting Issues 
Please report any issue regarding the rICESat2Veg package to Dr. Silva (c.silva@ufl.edu)

# Citing rICESat2Veg
Silva,C.A; Hamamura,C.rICESat2Veg: An R Package for NASA's Ice, Cloud, and Elevation Satellite (ICESat-2) Data Processing and Visualization for Terrestrial Applications.version 0.0.1, accessed on November. 22 2023, available at: <https://CRAN.R-project.org/package=rICESat2Veg>

# Disclaimer
**rICESat2Veg package has not been developted by the GEDI team. It comes with no guarantee, expressed or implied, and the authors hold no responsibility for its use or reliability of its outputs.**

