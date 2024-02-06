#' @include class.icesat2.h5_ds_local.R
#' @import R6 reticulate hdf5r
#' @export
ICESat2.h5_local <- R6::R6Class("ICESat2.h5_local", list(
    inherit = "ICESat2.h5",
    h5 = NULL,
    initialize = function(h5) {
        if (inherits(h5, "character")) {
            self$h5 <- H5File$new(h5)
        } else {
            self$h5 <- h5
        }
    },
    ls = function() {
        self$h5$ls()$name
    },
    exists = function(path) {
        self$h5$exists(path)
    },
    close_all = function() {
        self$h5$close_all()
    }
))

"[[.ICESat2.h5_local" <- function(x, i = NULL) {
    res <- x$h5[[i]]
    if (inherits(res, "H5D")) {
        return(ICESat2.h5_ds_local$new(ds = res))
    } else {
        return(ICESat2.h5_local$new(h5 = res))
    }
}
