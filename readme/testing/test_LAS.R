
install.packages("lidR")
require(lidR)

require(lidR)
atl08_path<-"C:\\Users\\c.silva\\Documents\\rICESat2Veg\\inst\\exdata\\ATL08_20220401221822_01501506_005_01.h5"
atl08_h5<-ATL08read(atl08_path=atl08_path)
atl03_path<-"C:\\Users\\c.silva\\Documents\\rICESat2Veg\\inst\\exdata\\ATL03_20220401221822_01501506_005_01.h5"
atl03_h5<-ATL03read(atl03_path=atl03_path)

## join
atl03_atl08_dt<-ATL03_ATL08_join_dt(atl03_h5,atl08_h5, beam = "gt1l")

atl03_atl08_dt@dt<-na.omit(atl03_atl08_dt@dt)
head(atl03_atl08_dt@dt)


dt2<-atl03_atl08_dt@dt[,c("lon_ph","lat_ph","ph_h")]
#class(dt2)<-"data.table"

names(dt2)<-c("X","Y","Z")
dt2$Classification <-atl03_atl08_dt@dt$classed_pc_flag

require(terra)
v <- vect(dt2, c("X", "Y"), crs="+proj=longlat")
u <- terra::project(v, "+proj=utm +zone=23")
utm <- crds(u)
nrow(utm)
nrow(dt2)

dt2$X<-utm[,1]
dt2$Y<-utm[,2]


header = LASheader(dt2)

# Record an EPSG code
#epsg(header) <- 32618

las <- LAS(dt2, header)

las@header@PHB[["X scale factor"]] <- 0.01
las@header@PHB[["Y scale factor"]] <- 0.01
las@header@PHB[["Z scale factor"]] <- 0.01

#st_crs(las)<-"EPSG:32723"

st_crs(las)$Name


#tlas <- sf::st_transform(las, sf::st_crs(26918))

head(las@data)

writeLAS(las,"C:\\Users\\c.silva\\Documents\\rICESat2Veg\\lasICE.laz")


plot(las[1:1000,], color="Classification")