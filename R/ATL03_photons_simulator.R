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
#'                   package="ICESat2VegR")
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
ATL08_read <-function(atl08_path) {

}
