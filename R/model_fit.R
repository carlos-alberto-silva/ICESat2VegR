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


#' Create ee random forest regression model
#'
#' @param x [`data.frame-class`]
#'
#' @include model_tools.R
#' @export
model_fit <- function(x, y, method = "randomForest") {
  fc <- build_fc(x, y)
  result <- randomForestRegression(fc, as.list(names(y)), as.list(names(x)), nTrees = 100, nodesize = 1)
  prepend_class(result, "ee.Classifier")
  return(result)
}
