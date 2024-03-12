#' @useDynLib ICESat2VegR

#' @export
earthaccess <- NULL

#' @export
ee <- NULL


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

.onAttach <- function(lib, pkg){
  info <- utils::packageDescription("ICESat2VegR")
  if (is.null(info$Date)){info$Date="2024-03-06 UTC"}
  base::packageStartupMessage(
    paste('\n##----------------------------------------------------------------##\n',
          'ICESat2VegR package, version ', info$Version, ', Released ', info$Date, '\n',
          'This package is based upon work supported by the NASA ICESat-2 ',
          'under grants No. ****. \n',
          '##----------------------------------------------------------------##',
          sep="")
  )
}