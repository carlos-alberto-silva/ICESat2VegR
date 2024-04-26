#' The class representing the h5 file opened from local files.
#'
#' @details
#' The variants `_cloud` and `_local` allows all the other functions
#' to use generic calls using the same interface, with each class
#' implementation is provided accordingly.
#'
#' The regular usage does not require the user to work with those classes
#' as most other provided functions will actually give access to the most
#' common necessities for working with ICESat-2 data.
#'
#' @field h5 A pointer to [`hdf5r::H5File-class`] in case
#' the user wants to access features not implemented yet
#' @field beams The [`character-class`] vector of beams available for the granule.
#' @field strong_beams The [`character-class`] vector of strong beams calculated using orbit_info
#' @field weak_beams The [`character-class`] vector of weak beams calculated using orbit_info
#' @field isOpen A flag to indicate if the file pointer has already been closed
#
#'
#' @include class.icesat2.R class.icesat2.h5ds_local.R
#' @import R6 reticulate hdf5r
ICESat2.h5_local <- R6::R6Class("ICESat2.h5_local", list(
  h5 = NULL,
  beams = NULL,
  strong_beams = NULL,
  weak_beams = NULL,
  isOpen = TRUE,
  #' @description Direct initialization should not be used, it is handled by `ATL0X_read()`
  #'
  #' @param h5 the result of the [`ATLAS_dataFinder()`]
  #'
  #' @return The class object
  initialize = function(h5) {
    if (inherits(h5, "character")) {
      self$h5 <- H5File$new(h5, mode = "r")
      groups <- self$ls()
      self$beams <- grep("gt[1-3][lr]", groups, value = TRUE)
      separated_beams <- list(
        grep("l$", self$beams, value = TRUE),
        grep("r$", self$beams, value = TRUE)
      )

      if (self$exists("orbit_info/sc_orient")) {
        sc_orient <- self[["orbit_info/sc_orient"]][]
        if (sc_orient == 2) {
          warning("Cannot determine the strong and weak beams from sc_orient == 2")
        }
        self$weak_beams <- separated_beams[[sc_orient + 1]]
        self$strong_beams <- setdiff(self$beams, self$weak_beams)
      } else {
        warning("Can't determine strong and weak beams, no [['orbit_info/sc_orient']] information!")
      }
    } else {
      self$h5 <- h5
    }
    prepend_class(self, "icesat2.h5")
  },
  #' @description Lists the groups and datasets that are within current group
  #'
  #' @return List the groups and datasets within the current path
  ls = function() {
    self$h5$ls()$name
  },
  #' @description Lists all grouops recursively
  #'
  #' @param recursive [`logical-class`], default FALSE. If TRUE it will list
  #' groups recursively and return the full path.
  #'
  #' @return The [`character-class`] representing
  ls_groups = function(recursive = FALSE) {
    data.table::as.data.table(
      self$h5$ls(recursive = recursive)
    )[obj_type == "H5I_GROUP", name]
  },
  #' @description Lists the available attributes
  #'
  #' @return [`character-class`] vector of attributes available
  ls_attrs = function() {
    hdf5r::list.attributes(self$h5)
  },
  #' @description Get datasets as data.table with columns (name, dataset.dims, rank)
  #'
  #' @param recursive [`logical-class`], default FALSE. If TRUE recursively searches
  #' and returns the full path.
  #'
  #' @return A [`data.table::data.table`] with the columns (name, dataset.dims, rank)
  dt_datasets = function(recursive = FALSE) {
    data.table::as.data.table(
      self$h5$ls(recursive)
    )[obj_type == "H5I_DATASET", list(name, dataset.dims, dataset.rank)]
  },
  #' @description Checks if a supplied group/dataset exist within the H5
  #'
  #' @param path [`character-class`] with the path to test
  #'
  #' @return [`logical-class`] TRUE or FALSE if it exists or not
  exists = function(path) {
    self$h5$exists(path)
  },
  #' @description Read an attribute from h5
  #'
  #' @param attribute [`character-class`] the address of the attribute to open
  attr = function(attribute) {
    hdf5r::h5attr(self$h5, attribute)
  },
  #' @description Safely closes the h5 file pointer
  #'
  #' @param silent [`logical-class`], default TRUE. Will cast warning messages if silent = FALSE.
  #'
  #' @return Nothing, just closes the file
  close_all = function(silent = TRUE) {
    if (self$isOpen) {
      self$isOpen <- FALSE
      tryCatch(
        {
          self$h5$close_all()
        },
        error = function(err) {
          self$h5$close()
        },
        silent = TRUE
      )
    }
  },
  #' @description Prints the data in a friendly manner
  #'
  #' @param ... Added for compatibility purposes with generic signature
  #'
  #' @return Outputs information about object
  print = function(...) {
    cat("Class: ICESat2.h5_local\n")
    cat("=============================\n")
    cat(paste0(self$ls(), collapse = "\n"))
    cat("\n")
  }
))

#' @export
"[[.ICESat2.h5_local" <- function(x, i) {
  res <- x$h5[[i]]
  if (inherits(res, "H5D")) {
    return(ICESat2.h5ds_local$new(ds = res))
  }
  return(ICESat2.h5_local$new(h5 = res))
}
