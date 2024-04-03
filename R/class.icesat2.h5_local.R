#' @include class.icesat2.h5ds_local.R class-ICESat2.h5.R
#' @import R6 reticulate hdf5r
#' @export
ICESat2.h5_local <- R6::R6Class("ICESat2.h5_local", list(
  inherit = "ICESat2.h5",
  h5 = NULL,
  beams = NULL,
  isOpen = TRUE,
  initialize = function(h5) {
    if (inherits(h5, "character")) {
      self$h5 <- H5File$new(h5, mode = "r")
      groups <- self$ls()
      if (self$attr("short_name") == "ATL03") {
        beams <- grep("gt[1-3][lr]", groups, value = TRUE)
        beamsList <- list()
        for (beam in beams) {
          if ("geolocation" %in% self[[beam]]$ls()) {
            beamsList[[""]] <- beam
          }
        }
        self$beams <- as.character(beamsList)
      } else {
        self$beams <- grep("gt[1-3][lr]", groups, value = TRUE)
      }
    } else {
      self$h5 <- h5
    }
  },
  ls = function() {
    self$h5$ls()$name
  },
  ls_groups = function(recursive = FALSE) {
    data.table::as.data.table(
      self$h5$ls(recursive = recursive)
    )[obj_type == "H5I_GROUP", name]
  },
  ls_attrs = function() {
    hdf5r::list.attributes(self$h5)
  },
  dt_datasets = function(recursive = FALSE) {
    data.table::as.data.table(
      self$h5$ls(recursive)
    )[obj_type == "H5I_DATASET", list(name, dataset.dims, dataset.rank)]
  },
  exists = function(path) {
    self$h5$exists(path)
  },
  attr = function(attribute) {
    hdf5r::h5attr(self$h5, attribute)
  },
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
