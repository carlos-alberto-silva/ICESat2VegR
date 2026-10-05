# Try logging in on earthaccess
#
# @param persist Logical. If TRUE, it will persist the login credentials in .netrc file.
#
# @return Nothing, just try to login in earthaccess
#
#' @keywords internal
earthaccess_login <- function(persist = TRUE) {
  # Test if earthaccess was loaded
  if (is.null(earthaccess)) {
    tryCatch(
      {
        earthaccess <- reticulate::import("earthaccess")
      },
      error = function(e) {
        stop("Earth access could not be loaded from reticulate, 
please run ICESat2Veg_configure().")
      }
    )
  }


  configured_netrc <- Sys.getenv("NETRC", unset = "")
  has_netrc <- nzchar(configured_netrc) && file.exists(configured_netrc)
  has_environment_credentials <- nzchar(Sys.getenv("EARTHDATA_TOKEN", unset = "")) ||
    (nzchar(Sys.getenv("EARTHDATA_USERNAME", unset = "")) &&
      nzchar(Sys.getenv("EARTHDATA_PASSWORD", unset = "")))

  if (has_environment_credentials) {
    auth <- earthaccess$login(strategy = "environment")
  } else {
    if (!has_netrc) {
      configured_netrc <- earthdata_login()
      has_netrc <- file.exists(configured_netrc)
    }
    auth <- if (has_netrc) earthaccess$login(strategy = "netrc") else NULL
  }

  authenticated <- !is.null(auth) && isTRUE(reticulate::py_to_r(auth$authenticated))
  if (!authenticated && interactive()) {
    auth <- earthaccess$login(strategy = "interactive", persist = persist)
    authenticated <- isTRUE(reticulate::py_to_r(auth$authenticated))
  }
  if (!authenticated) {
    stop("Could not authenticate in NASA earthaccess,
    please verify your credentials.")
  }

  message("Successfully logged in NASA earthaccess.")
}
