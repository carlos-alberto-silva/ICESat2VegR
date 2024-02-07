h5dtypes <- list()
h5dtypes[["int8"]] <- hdf5r::h5types$H5T_STD_I8LE
h5dtypes[["int32"]] <- hdf5r::h5types$H5T_STD_I32LE
h5dtypes[["int16"]] <- hdf5r::h5types$H5T_STD_I16LE
h5dtypes[["int64"]] <- hdf5r::h5types$H5T_STD_I64LE
h5dtypes[["uint16"]] <- hdf5r::h5types$H5T_STD_U16LE
h5dtypes[["float32"]] <- hdf5r::h5types$H5T_IEEE_F32LE
h5dtypes[["float64"]] <- hdf5r::h5types$H5T_IEEE_F64LE

#' @import R6
#' @export
ICESat2.h5ds_cloud <- R6::R6Class("ICESat2.h5ds_cloud", list(
    ds = NULL,
    dims = NULL,
    chunk_dims = NULL,
    initialize = function(ds) {
        self$ds <- ds
        self$dims <- rev(as.numeric(self$ds$shape))
        self$chunk_dims <- rev(as.numeric(self$ds$chunks))
    },
    get_type = function() {
        h5dtypes[[as.character(self$ds$dtype)]]
    },
    get_fill_value = function() {
        self$ds$fillvalue
    },
    get_creation_property_list = function() {
        py_plist <- self$ds$id$get_create_plist()
        number_filters <- py_plist$get_nfilters()
        
        plist <- hdf5r::H5P_DATASET_CREATE$new()
        for (ii in seq_len(number_filters - 1)) {
            py_filter <- py_plist$get_filter(ii)
            filter_id <- py_filter[[1]]
            flags <- py_filter[[2]]
            client_data <- as.integer(py_filter[[3]])


            plist$set_filter(filter_id, flags, client_data)
        }
        return(plist)
    }
))

#' @export
"[.ICESat2.h5ds_cloud" <- function(x, i, ...) {
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
