INT_MAX <- 2147483647


# Helper function to build a ee Feature Collection based on x and y vectors
#
# NOT EXPORTED
build_fc <- function(x, y = NULL) {
  all_data <- cbind(x, y)
  ee$FeatureCollection(lapply(seq_len(nrow(x)), function(ii) ee$Feature(NULL, as.list(all_data[ii]))))
}

# Helper function to convert the ee Feature Collection back to a data.table.
# This is used after extracting the data from the rasters.
#
# NOT EXPORTED
ee_to_dt <- function(sampled) {
  # Get size safely
  idx <- `system:index` <- NA
  sz <- tryCatch(sampled$size()$getInfo(), error = function(e) 0L)
  if (is.null(sz) || sz == 0L) return(data.table::data.table())

  # Bring features to R without using ee module
  feats <- sampled$toList(sz)$getInfo()  # list of ee.Feature dicts
  props <- lapply(feats, function(f) f$properties)

  # Bind with fill to handle missing props across features
  dt <- data.table::rbindlist(props, use.names = TRUE, fill = TRUE)

  # Build idx in R from system:index (no ee ops)
  if ("system:index" %in% names(dt)) {
    # Example: "123_0" -> 123 + 1
    dt[, idx := as.integer(gsub("_0.*", "", `system:index`)) + 1L ]
    dt[, `system:index` := NULL]
  }

  dt
}

# ee_to_dt <- function(sampled) {
#   sampled <- sampled$map(function(x) {
#     id <- x$getString("system:index")
#     id <- ee$Number(id$replace("_0", "")$decodeJSON())
#     x$set(list(
#       "idx" = id$add(1)
#     ))
#   })
#   columns <- sampled$first()$propertyNames()
#   columns <- columns$remove("system:index")
#
#   all_list <- list()
#   sampled$size()$getInfo()
#   nested_list <- sampled$reduceColumns(ee$Reducer$toList(columns$size()), columns)$values()$get(0)
#   dt <- data.table::rbindlist(nested_list$getInfo())
#   names(dt) <- columns$getInfo()
#
#   return(dt)
# }



#' @keywords internal
.get_pkg_module <- function() {
  Rcpp::Module("icesat2_module", PACKAGE = "ICESat2VegR")
}


#' Convert an R randomForest model to a Google Earth Engine randomForest classifier
#'
#' @description
#' Given a fitted [randomForest::randomForest] object, this function
#' serializes the forest using the internal Rcpp module `icesat2_module`
#' and constructs a corresponding Google Earth Engine random forest
#' **classifier** via
#' `ee$Classifier$decisionTreeEnsemble(rf_strings)`.
#'
#' This implementation always returns an `ee$Classifier` and does not use
#' `ee$Regressor`.
#'
#' @param rf A fitted [randomForest::randomForest] model object.
#'
#' @return
#' An Earth Engine `ee$Classifier` object created with
#' `ee$Classifier$decisionTreeEnsemble()` that can be used with
#' `ee$Image$classify()`.
#'
#' @export
build_ee_forest <- function(rf) {
  if (!inherits(rf, "randomForest")) {
    stop("'rf' must be a randomForest::randomForest object.")
  }

  # ---------------------------------------------------------------------------
  # 1) Get the Rcpp module and the buildForest function
  # ---------------------------------------------------------------------------
  mod <- tryCatch(
    Rcpp::Module("icesat2_module", PACKAGE = "ICESat2VegR"),
    error = function(e) {
      stop(
        "Could not load Rcpp module 'icesat2_module' from ICESat2VegR.\n",
        "Make sure the package is installed and its shared library is loaded.\n",
        "Original error: ", conditionMessage(e)
      )
    }
  )

  bf_fun <- tryCatch(
    mod$buildForest,
    error = function(e) {
      stop(
        "Rcpp module 'icesat2_module' does not expose 'buildForest'.\n",
        "Module summary:\n",
        paste(utils::capture.output(print(mod)), collapse = "\n")
      )
    }
  )

  # ---------------------------------------------------------------------------
  # 2) Serialize RF into a list/vector of tree strings
  # ---------------------------------------------------------------------------
  rf_strings <- tryCatch(
    bf_fun(rf),
    error = function(e) {
      stop("Error calling 'buildForest' in Rcpp module: ", conditionMessage(e))
    }
  )

  ok_list <- is.list(rf_strings) ||
    is.vector(rf_strings) ||
    inherits(rf_strings, "python.builtin.list")

  if (!ok_list) {
    stop(
      "buildForest(rf) must return a list/vector of serialized tree strings.\n",
      "Got object of class: ", paste(class(rf_strings), collapse = " ")
    )
  }

  # ---------------------------------------------------------------------------
  # 3) Build the Earth Engine classifier
  #    (assumes `ee` is already a valid Python Earth Engine module)
  # ---------------------------------------------------------------------------
  est <- tryCatch(
    ee$Classifier$decisionTreeEnsemble(rf_strings),
    error = function(e) {
      stop(
        "Failed to construct ee.Classifier.decisionTreeEnsemble from rf_strings.\n",
        "Check that your Earth Engine Python API supports this constructor.\n",
        "Original error: ", conditionMessage(e)
      )
    }
  )

  est
}
