#' @include class.icesat2.h5ds_local.R
#' @import R6 reticulate hdf5r
#' @export
ICESat2.h5_local <- R6::R6Class("ICESat2.h5_local", list(
    inherit = "ICESat2.h5",
    h5 = NULL,
    beams = NULL,
    initialize = function(h5) {
        if (inherits(h5, "character")) {
            self$h5 <- H5File$new(h5)
            groups <- self$ls()
            self$beams <- grep("gt[1-3][lr]", groups, value = TRUE)
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
    attr = function(attribute) {
        hdf5r::h5attr(self$h5, attribute)
    },
    close_all = function() {
        self$h5$close_all()
    }
))

#' @export
"[[.ICESat2.h5_local" <- function(x, i = NULL, j, ...) {
    res <- x$h5[[i]]
    if (inherits(res, "H5D")) {
        return(ICESat2.h5ds_local$new(ds = res))
    }
    return(ICESat2.h5_local$new(h5 = res))
}
