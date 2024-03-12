#' Given a stack image raster from GEE
#' retrieve the point geometry with values for the images
#'
#' @param stack A single image or a vector/list of images from Earth Engine.
#' @param geom A geometry from [`terra::SpatVector-class`] read with [`terra::vect`].
#' @param scale The scale in meters for the extraction (image resolution).
#'
#' @export
seg_gee_ancillary_dt_extract <- function(stack, geom, scale = 30) {
  sampled <- extract(stack, geom, scale)
  ee_to_dt(sampled)
}


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
  sampled <- sampled$map(function(x) {
    id <- x$getString("system:index")
    id <- ee$Number(id$replace("_0", "")$decodeJSON())
    x$set(list(
      "idx" = id$add(1)
    ))
  })
  columns <- sampled$first()$propertyNames()
  columns <- columns$remove("system:index")
  
  
  
  nested_list <- sampled$reduceColumns(ee$Reducer$toList(columns$size()), columns)$values()$get(0)
  dt <- data.table::rbindlist(nested_list$getInfo())
  names(dt) <- columns$getInfo()
  
  return(dt)
}

#' @export 
model_fit <- function(x, y, method = "randomForest") {
  fc <- build_fc(x, y)
  result <- randomForestRegression(fc, as.list(names(y)), as.list(names(x)), nTrees = 100, nodesize = 1)
  return(result)
}

#' @export
randomForestRegression <- function(
    featurecollection,
    property,
    train_properties,
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
      inputProperties = train_properties,
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
  if (!inherits(data, "ee.featurecollection.FeatureCollection")) {
    data <- build_fc(data, NULL)
  }
  
  predicted <- data$classify(x)
  classification <- predicted$aggregate_array("classification")$getInfo()
  return(classification)
}

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

build_fc <- function(x, y = NULL) {
  all_data <- cbind(x, y)
  ee$FeatureCollection(lapply(seq_len(nrow(x)), function(ii) ee$Feature(NULL, as.list(all_data[ii]))))
}

#' @export
var_select <- function(x, y, method = "forward", nboots = 10, nTrees = 100, train_split = 0.7, delta = 0.01) {
  ee_bandNames <- colnames(x)
  n <- nrow(x)
  train_size <- as.integer(round(n * train_split))
  validation_size <- n - train_size
  x <- build_fc(x, y)
  y <- names(y)

  selected_properties <- ee$List(list())
  rmses <- ee$List(list())

  minRmse <- ee_number(1e10)
  lastRmse <- ee_number(1e11)
  while (((1 - minRmse / lastRmse) >= delta)$getInfo() == 1) {
    result <- list(
      property = ee$List(list()),
      rmse = ee$List(list())
    )

    for (band in ee_bandNames) {
      current <- selected_properties$add(band)

      rmseList <- ee$List(list())
      message(gettextf("Testing %s", list(current$getInfo())), appendLF = TRUE)
      for (i in 1:nboots) {
        message(gettextf("\rBoot %d/%d", i, nboots), appendLF = FALSE)
        x <- x$randomColumn(seed = floor(runif(1) * INT_MAX))
        train_sample <- x$limit(train_size, "random")$select(current$add(y))
        validation_sample <- x$limit(validation_size, "random", ascending = FALSE)$select(current$add(y))


        randomForestClassifier <- randomForestRegression(train_sample, property = y, train_properties = current, nTrees = nTrees)
        classification <- validation_sample$classify(randomForestClassifier)
        classification2 <- classification$map(function(f) {
          f[["sqerror"]] <- (f[["h_canopy"]] - f[["classification"]]) ^ 2
          return(f)
        })
        rmse <- sqrt(classification2$aggregate_mean("sqerror"))
        rmseList[[]] <- rmse
      }
      message(appendLF = TRUE)

      meanRmse <- rmseList$reduce(ee$Reducer$mean())
      result$rmse <- result$rmse$add(meanRmse)
      result$property <- result$property$add(band)
    }
    lastRmse <- minRmse
    minRmse <- result$rmse$reduce(ee$Reducer$min())
    minIdx <- result$rmse$indexOf(minRmse)
    minRmse <- result$rmse$getNumber(minIdx)
    # ee$Dictionary(result)$getInfo()
    rmses <- rmses$add(minRmse)

    # rmses$getInfo()

    bandAdd <- result$property$getString(minIdx)
    selected_properties <- selected_properties$add(bandAdd)
    ee_bandNames <- ee_bandNames[ee_bandNames != bandAdd$getInfo()]
  }

  final_result <- ee$Dictionary(list(
    properties = selected_properties,
    rmse = rmses
  ))

  df <- data.frame(final_result$getInfo())
  return(df)
}
