#' The pointer to the `earthaccess` python reticulate module
#' @useDynLib ICESat2VegR
#' @import Rcpp Rdpack
#' @importFrom Rdpack reprompt
#' @export
earthaccess <- NULL

#' The pointer to the `earth-engine-api` python reticulate module
#'
#' @seealso https://developers.google.com/earth-engine/apidocs
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
        warning("Could not authenticate within earth-engine")
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
      "##  ICESat2VegR package, version ", info$Version, ", Released ", info$Date,
       "    #",
      "\n##----------------------------------------------------------------##\n",
      "\n",
      "This package is developed by Carlos A. Silva (c.silva@ufl.edu) and \n",
      "Caio Hamamura (caiohamamura@gmail.com) from the Forest Biometrics \n",
      "Remote Sensing and AI Lab (SilvaLab) at the University of Florida. \n",
      "ICESat2VegR is based upon work supported by NASA's Ice, Cloud, and \n",
      "Land Elevation Satellite (ICESat2), Carbon Monitoring System (CMS), \n",
      "and Commercial Smallsat Data Scientific Analysis (CSDSA) under \n",
      "grants No. 80NSSC23K0941, 80NSSC23K1257, and 80NSSC24K0055. \n",
      "\n",
      "Have a fantastic day filled with positivity and productivity! :) \n",
      "##----------------------------------------------------------------##",
      sep = ""
    )
  )
}
