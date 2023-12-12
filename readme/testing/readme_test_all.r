library(rICESat2Veg)
#atl08_path<-"C:\\Users\\c.silva\\Documents\\rICESat2Veg\\inst\\exdata\\ATL08_20220401221822_01501506_005_01.h5"
#atl08_h5<-ATL08_read(atl08_path=atl08_path)
#atl03_path<-"C:\\Users\\c.silva\\Documents\\rICESat2Veg\\inst\\exdata\\ATL03_20220401221822_01501506_005_01.h5"
#atl03_h5<-ATL03_read(atl03_path=atl03_path)

atl08_path<-"..\\inst\\exdata\\ATL08_20220401221822_01501506_005_01.h5"
atl03_path<-"..\\inst\\exdata\\ATL03_20220401221822_01501506_005_01.h5"

#atl08_path<-"C:\\Users\\c.silva\\Documents\\rICESat2Veg\\inst\\exdata\\ATL08_20220401221822_01501506_005_01.h5"
atl08_h5<-ATL08_read(atl08_path=atl08_path)
#atl03_path<-"C:\\Users\\c.silva\\Documents\\rICESat2Veg\\inst\\exdata\\ATL03_20220401221822_01501506_005_01.h5"
atl03_h5<-ATL03_read(atl03_path=atl03_path)


## join
atl03_atl08_dt<-ATL03_ATL08_join_dt(atl03_h5,atl08_h5, beam = "gt1l")

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
head(atl03_atl08_dt@dt)


# Clipping ATL08 Canopy Height Metrics by Geometry
# polygon_filepath <- "C:\\Users\\c.silva\\Documents\\rICESat2Veg\\inst\\exdata\\polygon.shp"
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

plot(atl03_atl08_dt@dt$lon_ph,atl03_atl08_dt@dt$lat_ph)
plot(polygon, add=T)

head(atl03_atl08_dt@dt)

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
    M1 = mean(x,na.rm = T),
    M2 = e1071::moment(x, order = 2, center = TRUE, na.rm = T) * length(x),
    M3 = e1071::moment(x, order = 3, center = TRUE, na.rm = T) * length(x),
    M4 = e1071::moment(x, order = 4, center = TRUE, na.rm = T) * length(x),
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
library(rICESat2Veg)
library(data.table)
#'
#'# Specifying the path to ATL08 file (zip file)
#'outdir = tempdir()
#'
#'atl08_zip <- system.file("extdata",
                  #"ATL08_20220401221822_01501506_005_01.zip",
                  #package="rICESat2Veg")
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
