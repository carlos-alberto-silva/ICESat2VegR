.datatable.aware = TRUE

.onUnload <- function(libpath) {
  library.dynam.unload("rICESat2Veg", libpath)
  invisible()
}

#' @useDynLib rICESat2Veg, .registration=TRUE
#' @export CPP_GDALDataset CPP_GDALRasterBand create_dataset
.onLoad <- function(libname, pkgname) {
  .RICESat2_CACHE <- new.env(FALSE, parent = globalenv())
  assign("old.PROJ_LIB", Sys.getenv("PROJ_LIB"), envir = .RICESat2_CACHE)
  Rcpp::loadModule("gdal_module", TRUE, TRUE)
  InitializeGDAL()
}

.onAttach <- function(libname, pkgname) {
  # Setup copied from rgdal package
  Sys.setenv("PROJ_LIB" = system.file("proj", package = pkgname)[1])
}
