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
        old_plist <- self$ds$get_create_plist()
        number_filters <- old_plist$get_nfilters()

        plist <- hdf5r::H5P_DATASET_CREATE$new()
        for (ii in (seq_len(number_filters) - 1)) {
            old_filter <- old_plist$get_filter(ii)
            filter_id <- old_filter$filter
            flags <- old_filter
            client_data <- as.integer(old_filter[[3]])

            plist$set_filter(filter_id, flags, client_data)
        }
        plist$set_fill_value(self$get_type(), self$get_fill_value())
        return(plist)
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
