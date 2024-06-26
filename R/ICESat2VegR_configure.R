check_conda <- function() {
  tryCatch(
    {
      reticulate::conda_version()
    },
    error = function(err) {}
  )
}

tryRestart <- function(auto = TRUE) {
  resp <- readline("I will need to restart RSession, would you like to proceed? y/[n]")
  if (substr(tolower(resp), 1, 1) == "y") {
    if (requireNamespace("rstudioapi", quietly = TRUE)) {
      try(rstudioapi::restartSession(), silent = TRUE)
    }
  }
  message("Please restart your R environment manually and run this function again.")
}

#' Configure environment for using Google Earth Engine functions
#' and accessing earthaccess cloud data
#'
#' @details
#' This function will configure the python environment as required by this package for:
#'  - Accessing the h5 files directly from the cloud for cloud computing
#'  - Using earth engine functions wrappers that are made available within this package.
#'
#' @return TRUE if everything works as expected, it is just an interactive installer
#' @export
ICESat2VegR_configure <- function() {
  if (is.null(check_conda())) {
    answer <- readline("Ananconda seems not to be available, would you like to install it? [y/N]: ")
    if (tolower(answer) == "y") {
      reticulate::install_miniconda()
      message("Miniconda was successfully installed!")
    }
  }

  if (!is.null(check_conda())) {
    reticulate::use_condaenv("r-reticulate")
    if (reticulate::py_module_available("h5py") == FALSE) {
      reticulate::py_install("h5py")
    }
    if (reticulate::py_module_available("h5py") == FALSE) {
      message("h5py is not compatible with this version of Python, upgrading...")
      reticulate::py_install("python")
      tryRestart()
    }

    if (reticulate::py_module_available("earthaccess") == FALSE) {
      reticulate::py_install("earthaccess")
    }
    if (reticulate::py_module_available("earthaccess") == FALSE) {
      message("earthaccess is not compatible with this version of Python, upgrading...")
      reticulate::py_install("python")
      tryRestart()
    }
    if (reticulate::py_module_available("ee") == FALSE) {
      reticulate::py_install("earthengine-api")
    }
  } else {
    warning("Something went wrong, conda is still not installed properly!")
  }

  if (!requireNamespace("ICESat2VegR", quietly = TRUE) || is.null(ee) || is.null(h5py) || is.null(earthaccess)) {
    warning("The package could not be loaded properly please read the readme in")
    warning("https://github.com/carlos-alberto-silva/ICESat2VegR and file a new issue")
    warning("if the problem is not solved.")
    return(FALSE)
  } else {
    return(TRUE)
  }
}

tryInitializeCloudCapabilities <- function() {
  # Reload namespace
  unloadNamespace("ICESat2VegR")
  loadNamespace("ICESat2VegR")

  if (is.null(h5py) || is.null(earthaccess)) {
    warning("The package's earthaccess functions were not loaded properly")
    warning("use ICESat2VegR_configure() function and reload the package!")
  }
}

tryInitializeEarthEngine <- function() {
  # Reload namespace
  unloadNamespace("ICESat2VegR")
  loadNamespace("ICESat2VegR")

  if (is.null(ee)) {
    warning("The package's earth-engine functions were not loaded properly")
    warning("use ICESat2VegR_configure() function and reload the package!")
  }

  tryCatch(
    {
      ee$Initialize()
    },
    error = function(e) {
      tryCatch(
        {
          ee$Authenticate()
        },
        error = function(e) {
          stop("Could not authenticate within earth-engine")
        }
      )
    }
  )
}
