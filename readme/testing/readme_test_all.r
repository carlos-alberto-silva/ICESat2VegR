# devtools::load_all()
library(ICESat2VegR)

# Specifying bounding box coordinates
lower_left_lon <- -96.0
lower_left_lat <- 40.0
upper_right_lon <- -100
upper_right_lat <- 42.0


# Specifying the date range
daterange <- c("2021-10-02", "2021-10-03")

# Extracting the path to ICESat-2 ATLAS data for the specified boundary box coordinates
atl08_granules_list_d <- ICESat2_dataFinder(
  short_name = "ATL08",
  lower_left_lon,
  lower_left_lat,
  upper_right_lon,
  upper_right_lat,
  version = "006",
  daterange = daterange,
  persist = TRUE,
  cloud_hosted = TRUE,
  cloud_computing = FALSE
)

atl08_granules_list_c <- ICESat2_dataFinder(
  short_name = "ATL08",
  lower_left_lon,
  lower_left_lat,
  upper_right_lon,
  upper_right_lat,
  version = "006",
  daterange = daterange,
  persist = TRUE,
  cloud_hosted = FALSE,
  cloud_computing = TRUE
)
ATL08_read(atl08_granules_list_c[1])[['orbit_info/sc_orient']][]

# Extracting the path to ICESat-2 ATLAS data for the specified boundary box coordinates
atl03_granules_list_d <- ICESat2_dataFinder(
  short_name = "ATL03",
  lower_left_lon,
  lower_left_lat,
  upper_right_lon,
  upper_right_lat,
  version = "006",
  daterange = daterange,
  persist = TRUE,
  cloud_hosted = TRUE,
  cloud_computing = FALSE
)

atl03_granules_list_c <- ICESat2_dataFinder(
  short_name = "ATL03",
  lower_left_lon,
  lower_left_lat,
  upper_right_lon,
  upper_right_lat,
  version = "006",
  daterange = daterange,
  persist = TRUE,
  cloud_hosted = FALSE,
  cloud_computing = TRUE
)


# reading ATL03 and ATL04 data from the cloud
atl03_cloud <- ATL03_read(atl03_granules_list_c[1])
atl08_cloud <- ATL08_read(atl08_granules_list_c[1])


# set working directory
output_dir <- tempdir()
setwd(output_dir)



atl08_path <- system.file("extdata/atl08_clip.h5", package = "ICESat2VegR")
atl03_path <- system.file("extdata/atl03_clip.h5", package = "ICESat2VegR")

# atl08_path<-"Z:\\01_Projects\\04_NASA_ICESat2\\11_others\\rICESat2Veg\\inst\\exdata\\ATL08_20220401221822_01501506_006_02.h5"
# atl03_path<-"Z:\\01_Projects\\04_NASA_ICESat2\\11_others\\rICESat2Veg\\inst\\exdata\\ATL03_20220401221822_01501506_006_02.h5"

# atl08_path<-"..\\inst\\exdata\\ATL08_20220401221822_01501506_005_01.h5"
# atl03_path<-"..\\inst\\exdata\\ATL03_20220401221822_01501506_005_01.h5"


## Reading ATL03 and ATL08 h5 files
```r
atl03_h5 <- ATL03_read(atl03_path = atl03_path)
atl08_h5 <- ATL08_read(atl08_path = atl08_path)
```

## Extracting ATL03 and ATL08 photon attributes
```r
atl03_photons_dt <- ATL03_photons_attributes_dt(atl03_h5 = atl03_h5, beam = "gt1r")
head(atl03_photons_dt)

##       lon_ph   lat_ph     h_ph quality_ph solar_elevation dist_ph_along
##        <num>    <num>    <num>      <num>           <num>         <num>
## 1: -103.7601 59.46874 325.9640          0        23.50526     0.3446315
## 2: -103.7601 59.46874 278.7379          0        23.50526     0.1905272
## 3: -103.7601 59.46873 429.6898          0        23.50526     1.3921306
## 4: -103.7601 59.46872 416.0034          0        23.50526     2.0568457
## 5: -103.7601 59.46872 405.5378          0        23.50526     2.0235429
## 6: -103.7601 59.46872 382.8574          0        23.50526     1.9488899

atl08_photons_dt <- ATL08_photons_attributes_dt(atl08_h5 = atl08_h5, beam = "gt1r")
head(atl08_photons_dt)

##    ph_segment_id   beam classed_pc_indx classed_pc_flag       ph_h d_flag delta_time
##            <int> <char>           <int>           <int>      <num>  <int>      <num>
## 1:        671214   gt1r               9               0 -1.2785034      1  134086702
## 2:        671214   gt1r              11               0 -1.0812683      1  134086702
## 3:        671214   gt1r              22               0 -2.2971802      1  134086702
## 4:        671214   gt1r              34               1 -0.3422546      1  134086702
## 5:        671214   gt1r              40               2  2.0543823      1  134086702
## 6:        671214   gt1r              61               2  0.9252319      1  134086702
```

## Extracting ATL08-derived terrain and canopy attributes by segments
```r
# Extracting
atl08_seg_att_dt <- ATL08_seg_attributes_dt(atl08_h5 = atl08_h5)
head(atl08_seg_att_dt)
class(atl08_seg_att_dt)

##    latitude longitude   beam h_canopy canopy_openness h_te_mean terrain_slope
##       <num>     <num> <char>    <num>           <num>     <num>         <num>
## 1: 59.46847 -103.7601   gt1r 3.259125       0.7998431  346.7855    0.02486196
## 2: 59.46758 -103.7603   gt1r 4.874969       1.1096631  346.1675   -0.06330263


summary(atl08_seg_att_dt)
##    latitude       longitude          beam              h_canopy     canopy_openness    h_te_mean     terrain_slope
## Min.   :59.47   Min.   :-103.8   Length:2           Min.   :3.259   Min.   :0.7998   Min.   :346.2   Min.   :-0.063303
## 1st Qu.:59.47   1st Qu.:-103.8   Class :character   1st Qu.:3.663   1st Qu.:0.8773   1st Qu.:346.3   1st Qu.:-0.041261
## Median :59.47   Median :-103.8   Mode  :character   Median :4.067   Median :0.9548   Median :346.5   Median :-0.019220
## Mean   :59.47   Mean   :-103.8                      Mean   :4.067   Mean   :0.9548   Mean   :346.5   Mean   :-0.019220
## 3rd Qu.:59.47   3rd Qu.:-103.8                      3rd Qu.:4.471   3rd Qu.:1.0322   3rd Qu.:346.6   3rd Qu.: 0.002821
## Max.   :59.47   Max.   :-103.8                      Max.   :4.875   Max.   :1.1097   Max.   :346.8   Max.   : 0.024862


# plotting histograms
atl08_seg_att_dt$h_canopy[atl08_seg_att_dt$h_canopy > 80] <- NA # set NA to values > 80 m
atl08_seg_att_dt$h_te_mean[atl08_seg_att_dt$h_te_mean > 5000] <- NA # set NA to values > 5000 m

oldpar <- par() # Save graphical parameters

par(mfrow = c(1, 2))
hist(atl08_seg_att_dt$h_canopy, col = "green", xlab = "height (m)", main = "h_canopy")
hist(atl08_seg_att_dt$h_te_mean, col = "#bd8421", xlab = "Elevation (m)", main = "h_te_mean")

par(oldpar) # Restore old graphical parameters

# Plotting ATL08 attribute on a map by segment coordinates
library(leaflet) # loading leaflet package
library(leafsync) # loading leafsync package

# Set breaks values for h_canopy and h_te_mean
options(scipen = 999)
h_canopy.breaks <- seq(0, 60, 10)
h_te_mean.breaks <- seq(3000, 4000, 200)

h_canopy_bins <- h_canopy.breaks[cut(atl08_seg_att_dt$h_canopy,
  breaks = h_canopy.breaks, labels = FALSE
)]

h_te_mean_bins <- h_te_mean.breaks[cut(atl08_seg_att_dt$h_te_mean,
  breaks = h_te_mean.breaks, labels = FALSE
)]

# set color palette
h_canopy.pal <- colorBin(palette = "BrBG", domain = h_canopy_bins)
h_te_mean.pal <- colorBin(palette = "BrBG", domain = h_te_mean_bins)

# plotting h_canopy
h_canopy.map <- leaflet() %>%
  setView(lng = -107.11, lat = 37.51314, zoom = 13.49) %>%
  addCircleMarkers(atl08_seg_att_dt$longitude,
    atl08_seg_att_dt$latitude,
    radius = 0.5,
    opacity = 1,
    color = h_canopy.pal(h_canopy_bins)
  ) %>%
  addScaleBar(options = list(imperial = FALSE)) %>%
  addProviderTiles(providers$Esri.WorldImagery) %>%
  addLegend(pal = h_canopy.pal, values = h_canopy_bins, title = "h_canopy (m)")

# plotting h_te_mean
h_te_mean.map <- leaflet() %>%
  setView(lng = -107.11, lat = 37.51314, zoom = 13.49) %>%
  addCircleMarkers(atl08_seg_att_dt$longitude,
    atl08_seg_att_dt$latitude,
    radius = 0.5,
    opacity = 1,
    color = h_te_mean.pal(h_te_mean_bins)
  ) %>%
  addScaleBar(options = list(imperial = FALSE)) %>%
  addProviderTiles(providers$Esri.WorldImagery) %>%
  addLegend(pal = h_te_mean.pal, values = h_te_mean_bins, title = "h_te_mean (m)")
sync(h_canopy.map, h_te_mean.map)
```

## Clipping ATL08 Terrain and Canopy Attributes

```r
# Clipping by  bounding box
# Define the bounding box
lower_left_lon=-107.7
lower_left_lat=42.75
upper_right_lon=-106.5
upper_right_lat=32.75

# Clip
atl08_seg_att_dt_clip <- ATL08_seg_attributes_dt_clipBox(atl08_seg_att_dt,
                                                         lower_left_lon,
                                                         upper_right_lon,
                                                         upper_right_lat,
                                                         lower_left_lat)

head(atl08_seg_att_dt_clip) # print the first six observations

##   latitude longitude beam  h_canopy canopy_openness h_te_mean  terrain_slope
## 1: 42.08012 -106.5000 gt1r  4.976807        1.069874  2075.281 -0.11553377658
## 2: 42.07922 -106.5001 gt1r  7.460938        1.273932  2063.163 -0.06796860695
## 3: 42.07832 -106.5003 gt1r 10.185303        2.310288  2056.079 -0.08571143448
## 4: 42.07742 -106.5004 gt1r  8.583984        2.444624  2050.611 -0.02976966091
## 5: 42.07652 -106.5005 gt1r  6.843018        2.086431  2048.529 -0.00003621415
## 6: 42.07563 -106.5006 gt1r  6.016113        1.481141  2049.141  0.00491552008

```

# Clipping by geometry

# Specify the path to shapefile
#poly_filepath <- system.file("extdata", "polygon.shp", package = "ICESat2VegR")
poly_filepath <- "Z:\\01_Projects\\04_NASA_ICESat2\\11_others\\rICESat2Veg\\inst\\exdata\\polygon.shp"

# Read shapefile
library(terra)
sppoly <- terra::vect(poly_filepath)
terra::crs(sppoly) = "epsg:4326"

class(atl08_seg_att_dt)

# Clip
atl08_seg_att_dt_clipg <- ATL08_seg_attributes_dt_clipGeometry(atl08_seg_att_dt, sppoly, split_by = "FID")
head(atl08_seg_att_dt_clipg) # print the first six observations

##   latitude longitude beam h_canopy canopy_openness h_te_mean terrain_slope nid  poly_id
## 1: 52.49692 -105.0109 gt1r 4.191498        1.060516  503.8196  0.0012361892   1       1
## 2: 52.47792 -105.0140 gt1r 7.399384        1.957333  513.7219  0.1044191048   2       1
## 3: 52.47721 -105.0141 gt1r 3.816345        1.160930  514.6264  0.0050510662   3       1
## 4: 52.46538 -105.0160 gt1r 6.513763        1.305076  504.8456 -0.0748032928   4       1
## 5: 52.46466 -105.0161 gt1r 6.364258        1.470707  503.5898  0.0159635153   5       1
## 6: 52.46377 -105.0163 gt1r 4.756866        1.288144  503.8241 -0.0004335443   6       1
```

## View ATL08 clipped data by bbox
m1<-leaflet() %>%
  setView(lng = -107.11, lat = 37.51314,zoom = 03) %>%
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
  addPolygons(data=sppoly,weight=1,col = 'white',
              opacity = 1, fillOpacity = 0) %>%
  addProviderTiles(providers$Esri.WorldImagery) %>%
  addLegend(pal = pal, values = atl08_seg_att_dt_clipg$poly_id,title ="Poly IDs" )
sync(m1, m2)

## Computing the top h_canopy at 30 m grid cell
atl08_seg_att_dt <- na.omit(atl08_seg_att_dt)
max_h_canopy <- ATL08_seg_attributes_dt_gridStat(atl08_seg_att_dt, func = max(h_canopy), res = 0.05)
plot(max_h_canopy,
  xlim = c(-107.2, -106.8), ylim = c(38, 39), col = viridis::inferno(8),
  main = "Max 'h_canopy'",
  xlab = "Langitude (degree)",
  ylab = "Latitude (degree)"
)

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

# Computing a series of h_canopy statistics at 30 m grid cellfrom customized function
atl08_seg_att_dt$h_canopy[atl08_seg_att_dt$h_canopy > 80] <- NA # set values > 80 m to NA m
h_canopy_metrics <- ATL08_seg_attributes_dt_gridStat(atl08_seg_att_dt, func = mySetOfMetrics(h_canopy), res = 0.05)
plot(h_canopy_metrics,
  xlim = c(-107.2, -106.8),
  ylim = c(38, 39),
  col = viridis::viridis(8),
  xlab = "Langitude (degree)",
  ylab = "Latitude (degree)"
)


# Define a list of metrics
multiple_metrics <- ATL08_seg_attributes_dt_gridStat(
  atl08_seg_att_dt,
  func = list(
    max_h_canopy = max(h_canopy),
    max_terrain_elevation = max(h_te_mean),
    mean_slope = mean(terrain_slope)
  ),
  res = 0.05
)
plot(multiple_metrics,
  xlim = c(-107.2, -106.8),
  ylim = c(38, 39),
  col = viridis::viridis(8),
  xlab = "Langitude (degree)",
  ylab = "Latitude (degree)"
)


## Merging ATL03 and ATL08 photon attributes
```r
atl03_atl08_dt<-ATL03_ATL08_photons_attributes_dt_join(atl03_h5, atl08_h5, beam = "gt1l")
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
```

## Plotting photon cloud
```r
# plot by "ph_h"
plot(atl03_atl08_dt, y="h_ph",colors= c("gray", "#bd8421", "forestgreen", "green"), beam = "gt1l",
     xlim = c(1500000, 1600000),pch=16)

# plot by "h_ph"
plot(atl03_atl08_dt, y="ph_h",colors= c("gray", "#bd8421", "forestgreen", "green"), beam = "gt1l",
     ylim=c(300,355),xlim=c(2400,3800), pch=16)
```

summary(atl03_atl08_dt[,2:3])

## Computing the mean of h_ph attribute at 0.05 degree grid cell
mean_h_ph <- ATL03_ATL08_photons_attributes_dt_gridStat(atl03_atl08_dt, func = mean(h_ph), res = 0.05)
plot(mean_h_ph,
     main="Mean h_ph",
     col=viridis::inferno(8),
     xlab="Langitude (degree)",
     ylab="Latitude (degree)")

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


## Clipping ATL03 and joined ATL03 and ATL08 photons attributes
```r
# Clipping by  bounding box

# ATL03 photons
atl03_photons_dt_clip <- ATL03_photons_attributes_dt_clipBox(atl03_photons_dt,
                                                             lower_left_lon,
                                                             upper_right_lon,
                                                             upper_right_lat,
                                                             lower_left_lat)
head(atl03_photons_dt_clip) # print the first six observations

# AT03 and ATL08 photons
atl03_atl08_photons_dt_clip <- ATL03_ATL08_photons_attributes_dt_clipBox(atl03_atl08_dt,
                                                                         lower_left_lon,
                                                                         upper_right_lon,
                                                                         upper_right_lat,
                                                                         lower_left_lat)
head(atl03_atl08_photons_dt_clip) # print the first six observations

# Clipping by geometry

# ATL03 photons
atl03_photons_dt_clipg <- ATL03_photons_attributes_dt_clipGeometry(atl03_photons_dt, sppoly, split_by = "FID")
head(atl03_photons_dt_clipg) # print the first six observations

# ATL03 and ATL08 photons
atl03_atl08_photons_dt_clipg <- ATL03_ATL08_photons_attributes_dt_clipGeometry(atl03_photons_dt, sppoly, split_by = "FID")
head(atl03_photons_dt_clipg) # print the first six observations
```

## Computing new segments a within length (e.g. 30m)
atl03_atl08_dt_seg <- ATL03_ATL08_segment_create(atl03_atl08_dt,
  segment_length = 30,
  centroid = "mean",
  output = NA,
  overwrite = FALSE
)

head(atl03_atl08_dt_seg)
## Computing terrain and canopy metrics by a within segment length
```r
# Computing the max canopy height at 30 m segments
max_ph_h_seg_30m <- ATL03_ATL08_compute_seg_attributes_dt_segStat(atl03_atl08_dt_seg, func=max(ph_h),
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

max_ph_h_seg_100m <-ATL03_ATL08_compute_seg_attributes_dt_segStat(atl03_atl08_dt, func=max(ph_h),
                                                                 ph_class=c(2,3),
                                                                 beam=c("gt1l", "gt1r", "gt2l", "gt2r", "gt3l", "gt3r"),
                                                                 quality_ph=0,
                                                                 night_flag=1)

head(max_ph_h_seg_100m)

##   segment max_ph_h latitude longitude dist_along
## 1:   1 2.025024 59.49107 -103.7541         15
## 2:   2 2.578125 59.49079 -103.7541         45
## 3:   3 3.153015 59.49051 -103.7542         75
## 4:   4 2.994751 59.49026 -103.7542        105
## 5:   5 2.130676 59.49000 -103.7543        135
## 6:   6 2.302063 59.48966 -103.7544        165

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
     sd = sd(x), # Sd of x
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

  ```r
 ## Statistics of ATL08 Terrain and Canopy Attributes by Geometry

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

























## Clipping joined ATL03 and ATL03 photons attributes

```r
# Clipping by  bounding box
```

## Clipping by geometry
```


# segment metrics
RH100max <-ATL03_ATL08_joined_dt_gridStat(atl03_atl08_dt, func=mean(ph_h),
                                res = 0.5,
                                ph_class=c(2,3),
                                beam=c("gt1l", "gt1r", "gt2l", "gt2r", "gt3l", "gt3r"),
                                quality_ph=0,
                                night_flag=1)

mySetOfMetrics <- function(x) {
  metrics <- list(
    min = min(x), # Min of x
    max = max(x), # Max of x
    mean = mean(x), # Mean of x
    sd = sd(x) # Sd of x
  )
  return(metrics)
}

canopy_metrics <-ATL08_canopy_dt_segStat(atl03_atl08_dt, func=mySetOfMetrics(ph_h),
                                seg_length = 30,
                                ph_class=c(2,3),
                                beam=c("gt1l", "gt1r", "gt2l", "gt2r", "gt3l", "gt3r"),
                                quality_ph=0,
                                night_flag=1)

Elev_max <-ATL08_terrain_dt_segStat(atl03_atl08_dt, func=mean(h_ph),
                                 seg_length = 100,
                                 ph_class=c(2,3),
                                 beam=c("gt1l", "gt1r", "gt2l", "gt2r", "gt3l", "gt3r"),
                                 quality_ph=0,
                                 night_flag=1)

windows()
plot(Elev_max@dt$latitude, Elev_max@dt$mean)

# extrac canopy attributes
atl08_canopy_dt<-ATL08_canopy_attributes_dt(atl08_h5=atl08_h5)
atl08_terrain_dt<-ATL08_terrain_attributes_dt(atl08_h5=atl08_h5)

head(atl08_canopy_dt@dt)
head(atl08_canopy_dt@dt)



# grid
ss<-ATL08_canopy_dt_gridStat(atl08_canopy_dt, func=max(h_canopy), res=0.5)
ss<-ATL08_canopy_dt_gridStat(atl08_canopy_dt, func=mySetOfMetrics(h_canopy), res=0.5)
ss
windows()
terra::plot(ss)
dev.off()

ss<-ATL08_terrain_dt_gridStat(atl08_terrain_dt, func=mySetOfMetrics(h_te_best_fit), res=0.5)
windows()
terra::plot(ss)
dev.off()

ss<-ATL03_ATL08_joined_dt_gridStat(atl03_atl08_dt, func=mySetOfMetrics(ph_h), res=0.5)
windows()
terra::plot(ss)
dev.off()
head(atl03_atl08_dt)


# Clipping ATL08 Canopy Height Metrics by Geometry
# polygon_filepath <- "C:\\Users\\c.silva\\Documents\\ICESat2VegR\\inst\\exdata\\polygon.shp"
polygon_filepath <- "..\\inst\\exdata\\polygon.shp"
library(terra)
polygon <- terra::vect(polygon_filepath)
polygon$FID<-c(1,2)

plot(polygon)

atl08_canopy_dt_clip <- ATL08_canopy_dt_clipGeometry(atl08_canopy_dt, polygon, split_by = "FID")
atl08_terrain_dt_clip <- ATL08_terrain_dt_clipGeometry(atl08_terrain_dt, polygon, split_by = "FID")
atl03_atl08_dt_clip <- ATL03_ATL08_joined_dt_clipGeometry(atl03_atl08_dt, polygon, split_by = "FID")

terrain_metrics_poly <-ATL08_terrain_dt_polyStat(atl08_terrain_dt_clip, func=mySetOfMetrics(h_te_best_fit),poly_id="poly_id")
h_canopy_metrics <-ATL08_canopy_dt_polyStat(atl08_canopy_dt_clip, func=mySetOfMetrics(h_canopy),id="poly_id")
max_ph_h <-ATL03_ATL08_joined_dt_polyStat(atl03_atl08_dt_clip, func=max(ph_h),poly_id="poly_id")


head(atl03_atl08_clip@dt)

head(atl03_atl08_clip)

plot(atl03_atl08_dt$lon_ph,atl03_atl08_dt$lat_ph)
plot(polygon, add=T)

head(atl03_atl08_dt)

plot(ss)

head(RH100max@dt)

plot(ss)


plot(ss)

head(canopy_metrics@dt)

ATL08_canopy_dt<-canopy_metrics

head(canopy_metrics@dt)

class(canopy_metrics)

head(canopy_metrics)


head(atl03_atl08_df)



#'
# Computing the maximum of RH100
RH100max <- polyStatsLevel2AM(level2AM_clip, func = max(rh100), id = NULL)
#'
max<-base::max



require(data.table)


head(atl03_atl08_df)






~data.table(
    n = length(x),
    M1 = mean(x,na.rm = TRUE),
    M2 = e1071::moment(x, order = 2, center = TRUE, na.rm = TRUE) * length(x),
    M3 = e1071::moment(x, order = 3, center = TRUE, na.rm = TRUE) * length(x),
    M4 = e1071::moment(x, order = 4, center = TRUE, na.rm = TRUE) * length(x),
    min = min(x, na.rm=T),
    max = max(x, na.rm=T)
  )

function(x1, x2) {
    combined = data.table()
    x1$n[is.na(x1$n)] = 0
    x1$M1[is.na(x1$M1)] = 0
    x1$M2[is.na(x1$M2)] = 0
    x1$M3[is.na(x1$M3)] = 0
    x1$M4[is.na(x1$M4)] = 0
    x1$max[is.na(x1$max)] = -Inf
    x1$min[is.na(x1$min)] = Inf

    combined$n = x1$n + x2$n
    delta = x2$M1 - x1$M1
    delta2 = delta * delta
    delta3 = delta * delta2
    delta4 = delta2 * delta2
#'
    combined$M1 = (x1$n * x1$M1 + x2$n * x2$M1) / combined$n
#'
    combined$M2 = x1$M2 + x2$M2 +
      delta2 * x1$n * x2$n / combined$n
#'
    combined$M3 = x1$M3 + x2$M3 +
      delta3 * x1$n * x2$n * (x1$n - x2$n) / (combined$n * combined$n)
    combined$M3 = combined$M3 + 3.0 * delta * (x1$n * x2$M2 - x2$n * x1$M2) / combined$n
#'
    combined$M4 = x1$M4 + x2$M4 + delta4 * x1$n * x2$n * (x1$n * x1$n - x1$n * x2$n + x2$n * x2$n) /
      (combined$n * combined$n * combined$n)
    combined$M4 = combined$M4 + 6.0 * delta2 * (x1$n * x1$n * x2$M2 + x2$n * x2$n * x1$M2) / (combined$n * combined$n) +
      4.0 * delta * (x1$n * x2$M3 - x2$n * x1$M3) / combined$n
#'
    combined$min = pmin(x1$min, x2$min, na.rm=F)
    combined$max = pmax(x1$max, x2$max, na.rm=F)
    return(combined)
}
list(
  sd = ~sqrt(M2/(n - 1)),
  skew = ~sqrt((n * (n - 1))) * ((sqrt(n) * M3) / (M2^1.5)) / (n - 2),
  kur = ~((n - 1) / ((n - 2) * (n - 3))) * ((n + 1) * ((n * M4) / (M2^2) - 3.0) + 6)
)

# Specifying the path to GEDI leveatl08_canopy_dt data (zip file)
library(ICESat2VegR)
library(data.table)
#'
#'# Specifying the path to ATL08 file (zip file)
#'outdir = tempdir()
#'
#'atl08_zip <- system.file("extdata",
                  #"ATL08_20220401221822_01501506_005_01.zip",
                  #package="ICESat2VegR")
#'
#'# Unzipping ATL08 file
#'atl08_path <- unzip(atl08_zip,exdir = outdir)
#'
# Bounding rectangle coordinates
ul_lat <- 59.50
ul_lon <- -108.3
lr_lat <- 26.99
lr_lon <- -103.8
res = 100 # meters
lat_to_met_factor = 1 / 110540
lon_to_met_factor = 1 / 111320
xres = lon_to_met_factor * res
yres = lat_to_met_factor * res
#'
agg_function = ~data.table(
    min = min(x),
    max = max(x),
    sum = sum(x),
    n = length(x))
#'
agg_join = function(agg1, agg2) {
agg1[is.na(agg1)] = 0
data.table(
    min = pmin(agg1$min, agg2$min),
    max = pmax(agg1$max, agg2$max),
    sum = agg1$sum + agg2$sum,
    n = agg1$n + agg2$n
)
}
#'
finalizer = list(
    mean = "sum/n",
    range = "max-min"
)
#'
outdir<-"..\\inst\\exdata\\"
list.files(outdir)
require(data.table)
ATL08_canopy_h5_gridStat(
  atl08_path = outdir,
  metrics = c("h_canopy"),
  out_root = file.path(outdir, "output"),
  lr_lat = lr_lat,
  ul_lat = ul_lat,
  lr_lon = lr_lon,
  ul_lon = ul_lon,
  res = c(xres, -yres),
  creation_options = c("COMPRESS=DEFLATE" ,
    "BIGTIFF=IF_SAFER",
    "TILED=YES",
    "BLOCKXSIZE=512",
    "BLOCKYSIZE=512"),
  agg_function = agg_function,
  agg_join = agg_join,
  finalizer = finalizer
  )
