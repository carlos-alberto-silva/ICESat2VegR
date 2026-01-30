#' Class representing dataset opened from locally using hdf5r
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
#' @keywords internal
ICESat2.h5ds_local <- R6::R6Class("ICESat2.h5ds_local", list(
  ds = NULL,
  dims = NULL,
  chunk_dims = NULL,

  #' @description Constructor using a dataset pointer
  #'
  #' @param ds Dataset pointer
  initialize = function(ds) {
    self$ds <- ds
    self$dims <- self$ds$dims
    self$chunk_dims <- self$ds$chunk_dims
  },

  #' @description Get the type of data for this dataset
  get_type = function() {
    self$ds$get_type()
  },

  #' @description Get the fill_value used by the dataset
  get_fill_value = function() {
    self$ds$get_fill_value()
  },


  #' @description Get creation property list for the dataset (plist)
  get_creation_property_list = function() {
    old_plist <- self$ds$get_create_plist()
    number_filters <- old_plist$get_nfilters()

    plist <- hdf5r::H5P_DATASET_CREATE$new()
    for (ii in (seq_len(number_filters) - 1)) {
      old_filter <- old_plist$get_filter(ii)
      filter_id <- old_filter$filter
      flags <- old_filter$flags
      client_data <- as.integer(old_filter[[3]])

      plist$set_filter(filter_id, flags, client_data)
    }
    plist$set_fill_value(self$get_type(), self$get_fill_value())
    return(plist)
  }
))

#' @keywords internal
"[.ICESat2.h5ds_local" <- function(x, ...) {
  args <- eval(substitute(alist(...)))
  for (ii in seq_along(args)) {
    jj <- args[[ii]]
    if (!missing(jj)) {
      args[[ii]] <- eval(jj, parent.frame())
    }
  }
  do.call("[", c(list(x$ds), args))
}

#' @keywords internal
"length.ICESat2.h5ds_local" <- function(x) {
  exp(sum(log(x$dims)))
}
