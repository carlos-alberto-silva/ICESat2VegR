#' @import R6
#' @export
ICESat2.h5ds_local <- R6::R6Class("ICESat2.h5ds_local", list(
    ds = NULL,
    dims = NULL,
    chunk_dims = NULL,
    initialize = function(ds) {
        self$ds <- ds
        self$dims <- self$ds$dims
        self$chunk_dims <- self$ds$chunk_dims
    },
    get_type = function() {
        self$ds$get_type()
    },
    get_fill_value = function() {
        self$ds$get_fill_value()
    },
    get_creation_property_list = function() {
        self$ds$get_create_plist()
    }
))

#' @export
"[.ICESat2.h5ds_local" <- function(x, i, ...) {
    try(
        {
            return(x$ds[i])
        },
        silent = TRUE
    )
    x$ds[]
}
