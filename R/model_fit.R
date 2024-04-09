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
#' @param x [`data.frame-class`]. The input features [`data.frame-class`] to use in [`randomForest::randomForest()`].
#' @param y [`numeric-class`]. The response vector to use in [`randomForest::randomForest()`]/
#' @param ... Additional parameters to pass to [`randomForest::randomForest()`] model.
#' @return A [`randomForest::randomForest`] object
#' 
#' @include model_tools.R
#' @import randomForest
#' @export
model_fit <- function(x, y, ...) {
  result <- randomForest::randomForest(x, y, ...)
  return(result)
}
