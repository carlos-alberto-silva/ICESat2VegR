#' Creates a tree ensamble Google Earth Engine string
#' from R's randomForest model
#'
#' @param rf [`randomForest::randomForest`]. A randomForest model
#' created by fitting the model with [`randomForest::randomForest()`].
#'
#' @details
#' This function will create
#'
#' @export
build_forest <- Rcpp::cppFunction('
CharacterVector build_forest(List rf) {
  void visitNode(std::stringstream&, List, int, int, int, int);
  int ntree = rf["ntree"];
  CharacterVector output(ntree);

  for (int ii = 0; ii < ntree; ii++) {
    std::stringstream temp;
    temp << "1) root 9999 9999 ";
    visitNode(temp, rf, 1, 1, 1, ii);
    output(ii) = temp.str();
  }

  return output;
}

void visitNode(std::stringstream& temp, List rf, int idx, int previdx = 1, int depth = 1, int tree = 1) {
  List rfobj = rf["forest"];
  // Accessing columns from DataFrame
  IntegerMatrix leftDaughter =  rfobj["leftDaughter"];
  int left_idx = leftDaughter(idx - 1, tree);

  NumericMatrix importance = rf["importance"];
  CharacterVector importanceNames = rownames(importance);
  IntegerMatrix bestvar = rfobj["bestvar"];
  int var_idx = bestvar(idx - 1, tree);
  auto split_var = var_idx > 0 ? Rcpp::as<std::string>(importanceNames[var_idx - 1]) : "";

  NumericMatrix xbestsplit = rfobj["xbestsplit"];
  double split_point = xbestsplit(idx - 1, tree);

  IntegerMatrix rightDaughter = rfobj["rightDaughter"];
  int right_idx = rightDaughter(idx - 1, tree);

  NumericMatrix nodepred = rfobj["nodepred"];
  double node_prediction = nodepred(idx - 1, tree);


  // Check if its a leaf node
  if (left_idx == 0) {
    temp << node_prediction << " *";
    return;
  } else {
    temp << node_prediction;
  }

  // Rule for the current node
  temp << std::endl << std::string(depth, \' \') << previdx * 2 << ") " <<
  split_var << ((left_idx % 2) == 0 ? " <= " : " > ") <<
  split_point << " 0 0 ";

  // Recursively visit left_idx and right_idx nodes
  visitNode(temp, rf, left_idx, previdx * 2, depth + 1, tree);
  temp << std::endl << std::string(depth, \' \') << previdx * 2 + 1 <<  ") " <<
  split_var << ((right_idx % 2) == 0 ? " <= " : " > ") <<
  split_point << " 0 0 ";
  visitNode(temp, rf, right_idx, previdx * 2 + 1, depth + 1, tree);
  return;
}')

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
  rf_strings <- build_forest(rf)
  ee$Classifier$decisionTreeEnsemble(rf_strings)
}