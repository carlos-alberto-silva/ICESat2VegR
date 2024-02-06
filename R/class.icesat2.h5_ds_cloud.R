#' @export
ICESat2.h5_ds_cloud <- R6Class("ICESat2.h5_ds_cloud", list(
    ds = NULL,
    initialize = function(ds) {
        self$ds <- ds
    },
    dims = function() {
        as.numeric(self$ds$shape)
    }
))

"[.ICESat2.h5_ds_cloud" <- function(x, i, ...) {
    try(
        {
            if (is.numeric(i)) {
                return(x$ds[(i - 1)])
            } else {
                return(x$ds[i])
            }
        },
        silent = TRUE
    )
    x$ds[]
}