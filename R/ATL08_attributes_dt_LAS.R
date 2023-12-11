

require(lidR)
atl08_path<-"C:\\Users\\c.silva\\Documents\\rICESat2Veg\\inst\\exdata\\ATL08_20220401221822_01501506_005_01.h5"
atl08_h5<-ATL08read(atl08_path=atl08_path)
atl03_path<-"C:\\Users\\c.silva\\Documents\\rICESat2Veg\\inst\\exdata\\ATL03_20220401221822_01501506_005_01.h5"
atl03_h5<-ATL03read(atl03_path=atl03_path)

## join
atl03_atl08_dt<-ATL03_ATL08join(atl03_h5,atl08_h5, beam = "gt1l")
dt2<-dt[,c("lon_ph","lat_ph","h_ph")]
class(dt2)<-"data.table"

names(dt2)<-c("X","Y","Z")
dt2$Classification <-dt$classed_pc_flag


header = LASheader(dt2)

# Record an EPSG code
#epsg(header) <- 32618


las <- LAS(dt2, header)

las@header@PHB[["X scale factor"]] <- 0.01
las@header@PHB[["Y scale factor"]] <- 0.01
las@header@PHB[["Z scale factor"]] <- 0.01



writeLAS(las,"C:\\Users\\c.silva\\Documents\\rICESat2Veg\\lasICE.laz")

plot(las)

#'Read ICESat-2 ATL08 data
#'
#'@description This function reads the ICESat-2 Land and
#'Vegetation Along-Track Products (ATL08) as h5 file.
#'
#'
#'@usage ATL08read(atl08_path)
#'
#'@param atl08_path File path pointing to ICESat-2 ATL08 data. Data in HDF5 Hierarchical Data Format (.h5).
#'
#'@return Returns an S4 object of class [`icesat2.atl08-class`] containing ICESat-2 ATL08 data.
#'
#'@seealso \url{https://icesat-2.gsfc.nasa.gov/sites/default/files/page_files/ICESat2_ATL08_ATBD_r006.pdf}
#'
#'@examples
#'# Specifying the path to ICESat-2 ATL08 data (zip file)
#'outdir = tempdir()
#'atl08_fp_zip <- system.file("extdata",
#'                   "ATL0802_A_2019108080338_O01964_T05337_02_001_01_sub.zip",
#'                   package="rICESat2Veg")
#'
#'# Unzipping ICESat-2 ATL08 data
#'atl08path <- unzip(atl08_fp_zip,exdir = outdir)
#'
#'# Reading ICESat-2 ATL08 data (h5 file)
#'atl08<-ATL08read(atl08_path=atl08_path)
#'
#'close(atl08)
#'@import hdf5r
#'@export
readATL08 <-function(atl08_path) {
  if (!is.character(atl08_path) | !tools::file_ext(atl08_path) == "h5") {
    stop("atl08_path must be a path to a h5 file")
  }

  atl08_h5 <- hdf5r::H5File$new(atl08_path, mode = 'r')
  atl08<- new("icesat2.atl08", h5 = atl08_h5)
  return(atl08)
}

