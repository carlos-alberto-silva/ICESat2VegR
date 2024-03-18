#' @useDynLib ICESat2VegR

#' @export
earthaccess <- NULL

#' @export
ee <- NULL


pkg_module <- Rcpp::Module("icesat2_module")

ee_cache <- new.env(parent = emptyenv())
ee_cache$search <- NULL

.onLoad <- function(libname, pkgname) {
  tryCatch(
    {
      reticulate::use_condaenv("r-reticulate")
    },
    error = function(err) {}
  )

  earthaccess <<- tryCatch(
    {
      reticulate::import("earthaccess", delay_load = TRUE, convert = FALSE)
    },
    error = function(err) {}
  )

  ee <<- tryCatch(
    {
      reticulate::import("ee", delay_load = TRUE, convert = TRUE)
    },
    error = function(err) {}
  )

  tryCatch(
    {
      ee$Initialize()
    },
    error = function(e) {
      tryCatch({ee$Authenticate()}, error = function(e){
        message("Could not authenticate within earth-engine")
      })
    }
  )
  h5py <<- tryCatch(
    {
      reticulate::import("h5py", delay_load = TRUE)
    },
    error = function(err) {}
  )


  if (is.null(h5py) || is.null(ee) || is.null(earthaccess)) {
    warning("The package earth-engine functions were not loaded properly")
    warning("use ICESat2VegR_configure() function and reload the package!")
  }
}


.onUnload <- function(libpath) {
  library.dynam.unload("ICESat2VegR", libpath)
}

.onAttach <- function(lib, pkg) {
  info <- utils::packageDescription("ICESat2VegR")
  if (is.null(info$Date)) {
    info$Date <- "2024-03-06 UTC"
  }
  base::packageStartupMessage(
    paste("\n##----------------------------------------------------------------##\n",
      "ICESat2VegR package, version ", info$Version, ", Released ", info$Date, "\n",
      "This package is based upon work supported by the NASA ICESat-2 ",
      "under grants No. ****. \n",
      "##----------------------------------------------------------------##",
      sep = ""
    )
  )
}
