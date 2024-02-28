#' Given a geometry with point samples and images from Earth Engine
#' retrieve the point geometry with values for the images
#'
#' @param images A single image or a vector/list of images from Earth Engine.
#' @param geom A geometry from [`terra::SpatVector-class`] read with [`terra::vect`].
#' @param scale The scale in meters for the extraction (image resolution).
#'
#' @export
extract <- function(images, geom, scale) {
  final <- c(images)
  tempjson <- tempfile(fileext = ".geojson")
  terra::writeVector(geom, tempjson, filetype = "geojson")
  parsed <- jsonlite::parse_json(readLines(tempjson))
  geojson <- ee$FeatureCollection(parsed)

  sampled <- final$sampleRegions(
    collection = geojson,
    scale = scale
  )

  return(sampled)
}

#' @export
ee_to_dt <- function(sampled) {
  ids <- sampled$aggregate_array("system:index")$getInfo()
  ids <- as.integer(gsub("_0", "", ids)) + 1
  df <- data.table::data.table(ids = ids)
  df[["ids"]] <- ids

  column_df <- sampled$first()$propertyNames()$getInfo()[-1]
  nested_list <- sampled$reduceColumns(ee$Reducer$toList(length(column_df)), column_df)$values()$get(0)

  pymain <- reticulate::import_main()
  pymain$nested_list <- nested_list
  reticulate::py_run_string("import numpy as np")
  reticulate::py_run_string("result = np.array(nested_list.getInfo())")


  df2 <- data.table::data.table(pymain$result)
  names(df2) <- column_df

  return(cbind(df, df2))
}

#' @export
randomForestRegression <- function(
    featurecollection,
    property_name,
    nTrees = 500,
    mtry = NULL,
    nodesize = 5,
    sampsize = .632,
    seed = floor((runif(1) * 2147483647 * 2) - 2147483647)) {
  if (is.null(mtry)) {
    mtry <- max(1, floor(length(featurecollection$first()$propertyNames()$getInfo()) / 3))
    classifier <- ee$Classifier$smileRandomForest(
      numberOfTrees = nTrees,
      variablesPerSplit = mtry,
      minLeafPopulation = nodesize,
      bagFraction = sampsize,
      seed = seed
    )$
      train(
      features = featurecollection,
      classProperty = "h_canopy"
    )$
      setOutputMode("REGRESSION")

    return(classifier)
  }

  ee$Classifier$smileRandomForest(
    numberOfTrees = nTrees,
    variablesPerSplit = mtry,
    minLeafPopulation = nodesize,
    bagFraction = sampsize,
    seed = seed
  )$
    train(
    features = featurecollection,
    classProperty = "h_canopy"
  )$
    setOutputMode("REGRESSION")
}


setRefClass("ee.Classifier")
#' @export
"predict.ee.Classifier" <- function(x, data, ...) {
  predicted <- data$classify(x)
  classification <- predicted$aggregate_array("classification")$getInfo()
  return(classification)
}

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
