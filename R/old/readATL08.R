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
#'@return Returns an S4 object of class ["icesat2.atl08_dt"] containing ICESat-2 ATL08 data.
#'
#'@seealso \url{https://icesat-2.gsfc.nasa.gov/sites/default/files/page_files/ICESat2_ATL08_ATBD_r006.pdf}
#'
#'@examples
#'# Specifying the path to ICESat-2 ATL08 data (zip file)
#'outdir = tempdir()
#'ATL08_fp_zip <- system.file("extdata",
#'                   "ATL0802_A_2019108080338_O01964_T05337_02_001_01_sub.zip",
#'                   package="rICESat2Veg")
#'
#'# Unzipping ICESat-2 ATL08 data
#'atl08_path <- unzip(ATL08_fp_zip,exdir = outdir)
#'
#'# Reading ICESat-2 ATL08 data (h5 file)
#'ATL08<-ATL08read(atl08_path=atl08_path)
#'
#'close(ATL08)
#'@import hdf5r
#'@export
ATL08read <-function(atl08_path) {
  ATL08_h5 <- hdf5r::H5File$new(atl08_path, mode = 'r')
  ATL08<- new("icesat2.atl08", h5 = ATL08_h5)
  return(ATL08)
}
