get_h5dtype <- function(dtype) {
  switch(dtype,
    "int8"    = hdf5r::h5types$H5T_STD_I8LE,
    "int16"   = hdf5r::h5types$H5T_STD_I16LE,
    "int32"   = hdf5r::h5types$H5T_STD_I32LE,
    "int64"   = hdf5r::h5types$H5T_STD_I64LE,
    "uint16"  = hdf5r::h5types$H5T_STD_U16LE,
    "float32" = hdf5r::h5types$H5T_IEEE_F32LE,
    "float64" = hdf5r::h5types$H5T_IEEE_F64LE
  )
}

#' Class representing dataset opened from the cloud using h5py
#'
#' @field ds Dataset pointer
#' @field dims Dimensions of the dataset
#' @field chunk_dims Chunk dimensions for the dataset
#'
#' @details
#' This class should not be instanced by the user, instead it is
#' automatically handled when the user opens a dataset using `[[]]`
#' operator from a h5 file.
#'
#' @import R6
#' @export
ICESat2.h5ds_cloud <- R6::R6Class("ICESat2.h5ds_cloud", list(
  ds = NULL,
  dims = NULL,
  chunk_dims = NULL,

  #' @description Constructor using a dataset pointer
  #'
  #' @param ds Dataset pointer
  initialize = function(ds) {
    self$ds <- ds
    self$dims <- rev(as.numeric(self$ds$shape))
    self$chunk_dims <- rev(as.numeric(self$ds$chunks))
  },

  #' @description Get the type of data for this dataset
  get_type = function() {
    get_h5dtype(as.character(self$ds$dtype))
  },

  #' @description Get the fill_value used by the dataset
  get_fill_value = function() {
    self$ds$fillvalue
  },

  
  #' @description Get creation property list for the dataset (plist)
  get_creation_property_list = function() {
    py_plist <- self$ds$id$get_create_plist()
    number_filters <- py_plist$get_nfilters()

    plist <- hdf5r::H5P_DATASET_CREATE$new()
    for (ii in (seq_len(number_filters) - 1)) {
      py_filter <- py_plist$get_filter(ii)
      filter_id <- py_filter[[1]]
      flags <- py_filter[[2]]
      client_data <- as.integer(py_filter[[3]])


      plist$set_filter(filter_id, flags, client_data)
    }
    plist$set_fill_value(self$get_type(), self$get_fill_value())
    return(plist)
  }
))

#' @export
"[.ICESat2.h5ds_cloud" <- function(x, ...) {
  args <- eval(substitute(alist(...)))
  for (ii in seq_along(args)) {
    jj <- args[[ii]]
    if (!missing(jj)) {
      args[[ii]] <- eval(jj, parent.frame())
    }
  }
  do.call("[", c(list(x$ds), args))
}

#' @export
"length.ICESat2.h5ds_cloud" <- function(x) {
  exp(sum(log(x$dims)))
}