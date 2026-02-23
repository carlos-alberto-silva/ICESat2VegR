#' The class representing the h5 file opened from the cloud for cloud computing
#'
#' @field h5 A pointer to `h5py` in case the user wants to access features not
#' implemented yet
#' @field beams The [`character-class`] vector of beams available for the granule.
#' @field strong_beams The [`character-class`] vector of strong beams calculated using orbit_info
#' @field weak_beams The [`character-class`] vector of weak beams calculated using orbit_info
#'
#' @details
#' Besides representing h5 files, it is also used to represent groups
#' opened with `[[]]`.
#'
#' The variants `_cloud` and `_local` allows all the other functions
#' to use generic calls using the same interface, with each class
#' implementation is provided accordingly.
#'
#' The regular usage does not require the user to work with those classes
#' as most other provided functions will actually give access to the most
#' common necessities for working with ICESat-2 data.
#'
#' @include class.icesat2.R class.icesat2.h5ds_cloud.R earthaccess.R
#' @import R6 reticulate
#' @export
ICESat2.h5_cloud <- R6::R6Class("ICESat2.h5_cloud", list(
  h5 = NULL,
  beams = NULL,
  strong_beams = NULL,
  weak_beams = NULL,
  #' @description
  #' Direct initialization should not be used, it is handled by `ATL0X_read()`
  #'
  #' @param h5 the result of the [`ATLAS_dataFinder()`] with `cloud_computing = TRUE`
  #'
  #' @return The class object
  initialize = function(h5) {
    if (inherits(h5, "icesat2.granules_cloud")) {
      stop("For now the package only works with one granule at a time
try with only one granule [i].")
    }
    if (inherits(h5, "icesat2.granule_cloud")) {
      os <- reticulate::import("os")
      sys <- reticulate::import("sys")
      py_builtins <- reticulate::import_builtins()
      sys$stdout <- py_builtins$open(os$devnull, "w")
      sys$stderr <- py_builtins$open(os$devnull, "w")

      # Try to refresh login
      earthaccess_login()
      granules <- earthaccess$open(list(h5))

      sys$stdout <- sys["__stdout__"]
      sys$stderr <- sys["__stderr__"]

      h5file <- h5py$File(granules[0])
      self$h5 <- h5file
      groups <- self$ls()
      self$beams <- grep("gt[1-3][lr]", groups, value = TRUE)
      separated_beams <- list(
        grep("l$", self$beams, value = TRUE),
        grep("r$", self$beams, value = TRUE)
      )

      sc_orient <- self[["orbit_info/sc_orient"]][]
      if (sc_orient == 2) {
        warning("Cannot determine the strong and weak beams from sc_orient == 2")
      }
      strong_index <- sc_orient + 1
      weak_index <- 2 - sc_orient

      self$weak_beams <- separated_beams[[weak_index]]
      self$strong_beams <- separated_beams[[strong_index]]
    } else {
      self$h5 <- h5
    }
    prepend_class(self, "icesat2.h5")
  },
  #' @description Lists the groups and datasets that are within current group
  #'
  #' @return List the groups and datasets within the current path
  ls = function() {
    pybuiltins <- reticulate::import_builtins()
    pybuiltins$list(self$h5$keys())
  },
  #' @description Lists all grouops recursively
  #'
  #' @param recursive [`logical-class`], default FALSE. If TRUE it will list
  #' groups recursively and return the full path.
  #'
  #' @return The [`character-class`] representing
  ls_groups = function(recursive = FALSE) {
    if (recursive) {
      pymain <- reticulate::import_main()
      pymain$keys <- list()
      pymain$h5 <- self$h5

      reticulate::py_run_string("
      import h5py
      h5.visit(lambda key: keys.append(key) if isinstance(h5[key], h5py.Group) else None)")
      all_keys <- pymain$keys

      return(all_keys)
    } else {
      all_items <- self$h5$ls()
    }
  },
  #' @description Lists the available attributes
  #'
  #' @return [`character-class`] vector of attributes available
  ls_attrs = function() {
    pybuiltins <- reticulate::import_builtins()
    pybuiltins$list(self$h5$attrs$keys())
  },
  #' @description Get datasets as data.table with columns (name, dataset.dims, rank)
  #'
  #' @param recursive [`logical-class`], default FALSE. If TRUE recursively searches
  #' and returns the full path.
  #'
  #' @return A [`data.table::data.table`] with the columns (name, dataset.dims, rank)
  dt_datasets = function(recursive = FALSE) {
    pybuiltins <- reticulate::import_builtins(convert = FALSE)
    nonrecursive <- function(callable) {
      items_visit <- reticulate::py_to_r(pybuiltins$list(self$h5$items()))
      for (x in items_visit) {
        callable(x[[1]], x[[2]])
      }
    }
    visitation <- if (recursive) self$h5$visititems else nonrecursive

    pyres <- pybuiltins$dict(
      name = pybuiltins$list(),
      dataset.dims = pybuiltins$list(),
      dataset.rank = pybuiltins$list()
    )

    visitation(function(x, obj) {
      if (inherits(obj, "h5py._hl.dataset.Dataset")) {
        pyres["name"]$append(x)
        pyres["dataset.dims"]$append(pybuiltins$str(" x ")$join(pybuiltins$map(pybuiltins$str, obj$shape)))
        pyres["dataset.rank"]$append(obj$ndim)
      }
    })
    res <- reticulate::py_to_r(pyres)
    data.table::as.data.table(data.frame(res))
  },
  #' @description Checks if a supplied group/dataset exist within the H5
  #'
  #' @param path [`character-class`] with the path to test
  #'
  #' @return [`logical-class`] TRUE or FALSE if it exists or not
  exists = function(path) {
    pymain <- reticulate::import_main()
    reticulate::py_run_string("def funcin(a, b): return a in b")
    pymain$funcin(path, self$h5)
  },
  #' @description Read an attribute from h5
  #'
  #' @param attribute [`character-class`] the address of the attribute to open
  attr = function(attribute) {
    pymain <- reticulate::import_main()
    pymain$temp <- self$h5$attrs
    pymain$attribute <- attribute
    reticulate::py_run_string("res = temp[attribute].decode('utf8')")$res
  },
  #' @description Safely closes the h5 file pointer
  #'
  #' @return Nothing, just closes the file
  close_all = function() {
    self$h5 <- NULL
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
"[[.ICESat2.h5_cloud" <- function(x, i = NULL) {
  res <- x$h5[[i]]
  if (inherits(res, "h5py._hl.dataset.Dataset")) {
    return(ICESat2.h5ds_cloud$new(ds = res))
  } else {
    return(ICESat2.h5_cloud$new(h5 = res))
  }
}
