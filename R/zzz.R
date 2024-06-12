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

# Private h5py module
h5py <- NULL

# Private package environment
pkg_env <- environment()

# Module with the Rcpp ANN indexing functions
pkg_module <- Rcpp::Module("icesat2_module")

# Private cache for the Earth Engine search
ee_cache <- new.env(parent = emptyenv())
ee_cache$search <- NULL

.onLoad <- function(libname, pkgname) {
  if (reticulate::py_module_available("earthaccess")) {
    earthaccess <<- reticulate::import("earthaccess", convert = FALSE)
  }

  if (reticulate::py_module_available("ee")) {
    ee <<- reticulate::import("ee", convert = TRUE)
    tryCatch(
      {
        ee$Initialize()
      },
      error = function(e) {
        tryCatch(
          {
            ee$Authenticate()
          },
          error = function(e) {}
        )
      }
    )
  }

  if (reticulate::py_module_available("h5py")) {
    h5py <<- reticulate::import("h5py", convert = TRUE)
  }

  loadGdal(pkgname)
}


loadGdal <- function(pkg_name) {
  Rcpp::loadModule("gdal_module", TRUE, TRUE)
  proj_path <- ""
  proj_version <- GetProjVersion()

  major <- proj_version[1]
  minor <- proj_version[2]

  # Since PROJ 9.1, the data files are in PROJ_DATA instead of PROJ_LIB
  if (major > 9 ||
    (major == 9 && minor >= 1)) {
    proj_path <- Sys.getenv("PROJ_DATA")
  } else {
    proj_path <- Sys.getenv("PROJ_LIB")
  }

  if (proj_path == "") {
    proj_path <- system.file("proj", package = pkg_name)[1]
  }

  InitializeGDAL(proj_path)
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
