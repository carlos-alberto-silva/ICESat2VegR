#' @export
earthaccess <- NA

#' @export
ee <- NA

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
