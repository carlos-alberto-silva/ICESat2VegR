#' @export
earthaccess <- NA


.onLoad <- function(libname, pkgname) {  
    earthaccess <<- reticulate::import("earthaccess", delay_load = TRUE, convert = FALSE)
    h5py <<- reticulate::import("h5py", delay_load = TRUE)
}