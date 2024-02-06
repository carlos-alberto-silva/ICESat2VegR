#' @import R6
#' @export
ICESat2.h5_ds_local <- R6::R6Class("ICESat2.h5_ds_local", list(
    ds = NULL,
    dims = NULL,
    initialize = function(ds) {
        self$ds <- ds
        self$dims <- self$ds$dims
    }
))

"[.ICESat2.h5_ds_local" <- function(x, i, ...) {
    try(
        {
            return(x$ds[i])
        },
        silent = TRUE
    )
    x$ds[]
}
