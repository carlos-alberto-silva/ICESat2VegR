#' @include class.icesat2.h5ds_cloud.R class-ICESat2.h5.R
#' @import R6 reticulate
#' @export
ICESat2.h5_cloud <- R6::R6Class("ICESat2.h5_cloud", list(
  inherit = "ICESat2.h5",
  h5 = NULL,
  beams = NULL,
  strong_beams = NULL,
  weak_beams = NULL,
  initialize = function(h5) {
    if (inherits(h5, "icesat2.granules_cloud")) {
      stop("For now the package only works with one granule at a time
try with only one granule [i].")
    }
    if (inherits(h5, "icesat2.granule_cloud")) {
      granules <- earthaccess$open(list(h5))
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
      self$weak_beams <- separated_beams[[sc_orient + 1]]
      self$strong_beams <- setdiff(self$beams, self$weak_beams)
    } else {
      self$h5 <- h5
    }
  },
  ls = function() {
    pybuiltins <- reticulate::import_builtins()
    pybuiltins$list(self$h5$keys())
  },
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
  exists = function(path) {
    pymain <- reticulate::import_main()
    pymain$h5 <- self$h5
    pymain$path <- path
    reticulate::py_run_string("res = path in h5")$res
  },
  attr = function(attribute) {
    pymain <- reticulate::import_main()
    pymain$temp <- self$h5$attrs
    pymain$attribute <- attribute
    reticulate::py_run_string("res = temp[attribute].decode('utf8')")$res
  },
  close_all = function() {
    self$h5 <- NULL
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
