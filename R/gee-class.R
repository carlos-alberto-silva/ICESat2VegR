library(R6)
eeCollection <- R6::R6Class(
  "eeCollection",
  public = list(
    collection_id = NA,
    bands = list(),
    expression = NA,
    collection = NA,
    initialize = function(collection_id, bands, expression) {
      self$collection_id <- collection_id
      self$bands <- bands
      if (exists("ee") && ee == reticulate::import("ee") && ee$Authenticate()) {
        self$collection <- ee$ImageCollection(collection_id)
        self$bands <- self$collection$first()$bandNames()$getInfo()
      }
    },
    count = function() {
      self$collection$size()$getInfo()
    },
    print = function(...) {
      cat("Bands: ")
      cat(self$bands, sep = " ")
      cat("\n")
      invisible(self)
    }
  ),
)

eeImage <- R6::R6Class(
  "eeImage",
  public = list(
    image = NA,
    bands = c(),
    initialize = function(img) {
      if (exists("ee") && ee == reticulate::import("ee") && ee$Authenticate()) {
        if (inherits(img, "character")) {
          self$image <- ee$Image(img)
        } else {
          self$image <- img
        }
        self$bands <- self$image$bandNames()$getInfo()
      }
    }
  )
)



"[[.eeCollection" <- function(x, ii) {
  eeImage$new(x$select(ii))
}
