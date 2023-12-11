#'Read ICESat-2 ATL03 data
#'
#'@description This function reads the ICESat-2 Global Geolocated Photons (ATL03) Product (ATL03) as h5 file.
#'
#'@usage readATL03(ATL03path)
#'
#'@param ATL03path File path pointing to ICESat-2 ATL03 data. Data in HDF5 Hierarchical Data Format (.h5).
#'
#'@return Returns an S4 object of class [`icesat2.ATL03-class`] containing ICESat-2 ATL03 data.
#'
#'@seealso \url{https://icesat-2.gsfc.nasa.gov/sites/default/files/page_files/ICESat2_ATL03_ATBD_r006.pdf}
#'
#'@examples
#'# Specifying the path to ICESat-2 ATL03 data (zip file)
#'outdir = tempdir()
#'ATL03_fp_zip <- system.file("extdata",
#'                   "ATL0302_A_2019108080338_O01964_T05337_02_001_01_sub.zip",
#'                   package="rICESat2Veg")
#'
#'# Unzipping ICESat-2 ATL03 data
#'ATL03path <- unzip(ATL03_fp_zip,exdir = outdir)
#'
#'# Reading ICESat-2 ATL03 data (h5 file)
#'ATL03<-readATL03(ATL03path=ATL03path)
#'
#'close(ATL03)
#'@import hdf5r
#'@export
readATL03 <-function(ATL03path) {
  ATL03_h5 <- hdf5r::H5File$new(ATL03path, mode = 'r')
  ATL03<- new("icesat2.ATL03", h5 = ATL03_h5)
  return(ATL03)
}

