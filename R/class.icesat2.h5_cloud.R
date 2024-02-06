#' @import R6 reticulate
#' @export
ICESat2.h5_cloud <- R6::R6Class("ICESat2.h5_cloud", list(
    inherit = "ICESat2.h5",
    h5 = NULL,
    initialize = function(h5) {
        self$h5 <- h5
    },
    ls = function() {
        pymain <- reticulate::import_main()
        pymain$temp <- self$h5$keys()
        reticulate::py_run_string("temp = list(temp)")$temp
    },
    exists = function(path) {
        pymain <- reticulate::import_main()
        pymain$h5 <- atl08_h5@h5
        pymain$path <- "gt1r"
        reticulate::py_run_string("res = path in h5")$res
    },
    close_all = function() {
        self$h5 <- NULL
    }
))

"[[.ICESat2.h5_cloud" <- function(x, i = NULL) {
    res <- x$h5[[i]]
    if (inherits(res, "h5py._hl.dataset.Dataset")) {
        return(ICESat2.h5_ds_cloud$new(ds = res))
    } else {
        return(ICESat2.h5_cloud$new(h5 = res))
    }
}
