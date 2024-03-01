#' @useDynLib ICESat2VegR

#' @export
earthaccess <- NA

#' @export
ee <- NA


pkg_module <- Rcpp::Module("icesat2_module")

ee_cache <- new.env(parent = emptyenv())
ee_cache$search <- NULL

.onLoad <- function(libname, pkgname) {
  earthaccess <<- reticulate::import("earthaccess", delay_load = TRUE, convert = FALSE)
  ee <<- reticulate::import("ee", delay_load = TRUE, convert = TRUE)
  tryCatch(
    {
      ee$Initialize()
    },
    error = function(e) {
      ee$Authenticate()
    }
  )
  h5py <<- reticulate::import("h5py", delay_load = TRUE)
}


.onUnload <- function (libpath) {
  library.dynam.unload("ICESat2VegR", libpath)
}