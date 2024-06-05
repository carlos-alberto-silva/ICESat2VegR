#' Try logging in on earthaccess
#'
#' @param persist Logical. If TRUE, it will persist the login credentials in .netrc file.
#'
#' @return Nothing, just try to login in earthaccess
#'
#' @examples
#' # Try to login in NASA earthaccess
#' \donttest{
#' # Shouldn't be tested because it relies on a web services requiring authentication.
#' earthaccess_login()
#' }
#'
#' @export
earthaccess_login <- function(persist = TRUE) {
  # Test if earthaccess was loaded
  if (is.null(earthaccess)) {
    tryCatch(
      {
        earthaccess <- reticulate::import("earthaccess")
      },
      error = function(e) {
        stop("Earth access could not be loaded from reticulate, please run ICESat2Veg_configure().")
      }
    )
  }

  auth <- earthaccess$login(strategy = "environment")
  if (!py_to_r(auth$authenticated) && file.exists(".netrc")) {
    auth <- earthaccess$login(strategy = "netrc")
  }
  if (!py_to_r(auth$authenticated)) {
    auth <- earthaccess$login(strategy = "interactive", persist = persist)
  }
  if (!py_to_r(auth$authenticated)) {
    stop("Could not authenticate in NASA earthaccess,
    please verify your credentials.")
  }

  if (py_to_r(auth$authenticated)) {
    message("Successfully logged in NASA earthaccess.")
  } else {
    stop("Could not authenticate in NASA earthaccess,
    please verify your credentials.")
  }
}
