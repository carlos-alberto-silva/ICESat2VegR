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
  sampled <- sampled$map(function(x) {
    id <- x$getString("system:index")
    id <- ee$Number(id$replace("_0", "")$decodeJSON())
    x$set(list(
      "idx" = id$add(1)
    ))
  })
  columns <- sampled$first()$propertyNames()
  columns <- columns$remove("system:index")

  all_list <- list()
  sampled$size()$getInfo()
  nested_list <- sampled$reduceColumns(ee$Reducer$toList(columns$size()), columns)$values()$get(0)
  dt <- data.table::rbindlist(nested_list$getInfo())
  names(dt) <- columns$getInfo()

  return(dt)
}


#' Given an R [`randomForest::randomForest()`] model, transform to a Google Earth Engine randomForest model
#' 
#' @param rf the [`randomForest::randomForest()`] model object.
#' 
#' @return The Google Earth Engine classifier
#' 
#' @export
build_ee_forest <- function(rf) {
  rf_strings <- pkg_module$buildForest(rf)
  ee$Classifier$decisionTreeEnsemble(rf_strings)
}