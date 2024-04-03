
devtools::install()
devtools::load_all()


 # Specifying bounding box coordinates
 lower_left_lon <- -96.0
 lower_left_lat <- 40.0
 upper_right_lon <- -100
 upper_right_lat <- 42.0

 # Specifying the date range
 daterange <- c("2019-07-01", "2020-05-22")

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
 atl03_cloud <- ATL03_read(atl03_granules_list_c[[1]])
 atl08_cloud <- ATL03_read(atl08_granules_list_c[[1]])

 # set working directory
 output_dir <- tempdir()
 setwd(output_dir)

 # download ATL03 and ATL04 data
 ICESat2_dataDownload(atl03_granules_list_d[1], ".")
 ICESat2_dataDownload(atl08_granules_list_d[1], ".")

 # List ATL03 and ATL08 data in the directory
 atl03_path <- list.files(pattern = "^ATL03")[1]
 atl08_path <- list.files(pattern = "^ATL08")[1]

 atl03_path<-"C:\\Users\\c.silva\\Documents\\ICESat2VegR\\readme\\testing\\processed_ATL03_20190711201711_02100406_005_01.h5"
 atl08_path<-"C:\\Users\\c.silva\\Documents\\ICESat2VegR\\readme\\testing\\processed_ATL08_20190711201711_02100406_005_01.h5"

  # Read ATL03 and ATL08 data
 atl03_h5 <- ATL03_read(atl03_path)
 atl08_h5 <- ATL08_read(atl08_path)


 atl03_h5[[paste0("gtr", "/geolocation/solar_elevation")]][]

 atl03_photons_dt<-ATL03_photons_attributes_dt(atl03_h5=atl03_h5, beam="gtr")
 head(atl03_photons_dt)


 atl03_atl08_dt <- ATL03_ATL08_photons_attributes_dt_join(atl03_h5, atl08_h5, beam = "gt1l")
 plot(atl03_atl08_dt, "h_ph", xlim = c(406000, 410000), ylim = c(180, 220), pch = 20, cex = 0.7)
