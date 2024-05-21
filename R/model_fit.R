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



#' Fit different models based on x and y inputs
#'
#' @param x input data.frame of predictors (independ variables) the data from
#' @param y input vector of observed data
#' @param method the model to use the data from. The options are "lm" for linear
#' regression, "rf" for random forest, "nnt" for neural network, "svm" for support
#' vector machine and knn.* for k-nearest neighbors (see [`yaImpute::yai()`] for * methods).
#' The default is "lm".
#' @param LOOCV [`logical-class`], default FALSE. Flag do determine if we should
#' compute leave-one-out cross validation.
#' @param ... Other parameters to pass to the model
#' @param size [`integer-class`], default 40, used for neural network number of
#' neurons only
#' @param linout [`logical-class`], default TRUE. Used only for neural network as
#' the activation function for the neuron to be linear instead of logistic.
#'
#' @return A list containing the method, model and cross validation data
#' if `LOOCV = TRUE`.
#' @export
model_fit <- function(x, y, method = "lm", LOOCV = FALSE, ..., size = 40, linout = TRUE) {
  oldpar <- graphics::par(no.readonly = TRUE)
  if (LOOCV) {
    on.exit(graphics::par(oldpar))
    graphics::par(mfrow = c(1, 2))
  }

  nnt_adapter <- function(y, x, ...) {
    if (requireNamespace("nnet", quietly = TRUE) == FALSE) {
      stop("Package nnet is not available for using nnt")
    }
    model <- nnet::nnet(y ~ ., data = x, ..., size = size, linout = linout)
    return(model)
  }

  yai_adapter <- function(y, x, ...) {
    if (requireNamespace("yaImpute", quietly = TRUE) == FALSE) {
      stop("Package yaImpute is not available for using yai")
    }
    yaImpute::yai(x = x, y = y, method = gsub("knn\\.", "", method), ...)
  }

  lm_adapter <- function(y, x, ...) {
    model <- stats::lm(formula = y ~ ., data = x, ...)
    return(model)
  }

  model_funs <- list(
    lm = lm_adapter,
    rf = tryCatch(
      randomForest::randomForest,
      error = function(x) stop("Package randomForest is not available")
    ),
    nnt = nnt_adapter,
    svm = tryCatch(
      e1071::svm,
      error = function(x) stop("Package e1071 is not available for using svm")
    )
  )

  model_fun <- yai_adapter

  if (method %in% names(model_funs)) {
    model_fun <- model_funs[[method]]
  }

  model <- model_fun(y = y, x = x, ...)

  results <- list()
  pred_model <- predict(model)
  pred_model <- if ("y" %in% names(pred_model)) pred_model[["y"]] else pred_model
  y_names <- if (is.null(colnames(y))) "" else colnames(y)
  if (length(dim(y)) == 2) {
    on.exit(graphics::par(oldpar))
    graphics::par(mfrow = c(dim(y)[2], graphics::par()$mfrow[2]))
  }
  for (yy in y_names) {
    y_local <- if (yy == "") y else y[[yy]]
    pred_model2 <- if (yy == "") pred_model else pred_model[[yy]]
    pred_model2 <- cbind(method = rep(method, nrow(x)), obs = y_local, pred = pred_model2)

    main_title <- if (yy == "") paste("model", method) else paste("model:", method, "var:", yy)
    stats_model <- cbind(
      method = rep(method, 6),
      stats_model(y_local, as.numeric(pred_model2[, 3]), main = main_title)
    )
    if (LOOCV == TRUE) {
      pb <- utils::txtProgressBar(min = 0, max = nrow(x), style = 3)
      predloocv <- NULL
      for (i in 1:nrow(x)) {
        utils::setTxtProgressBar(pb, i)
        model_i <- model_fun(y = y_local[-i], x = x[-i, ], ...)

        pred_i <- predict(model_i, x[i, ])

        predloocv <- rbind(predloocv, cbind(method = method, obs = y_local[i], pred = pred_i))
      }
      stats_loocv <- cbind(
        method = rep(method, 6),
        stats_model(as.numeric(predloocv[, 2]), as.numeric(predloocv[, 3]), main = paste("LOOCV", method))
      )
      colnames(pred_model2) <- c("method", "obs", "pred")
      colnames(predloocv) <- colnames(pred_model2)
      colnames(stats_model) <- c("method", "stat", "value")
      colnames(stats_loocv) <- c("method", "stat", "value")
      close(pb)
      results[[yy]] <- list(
        method = method, model = model, predModel = as.data.frame(pred_model2),
        predLOOCV = as.data.frame(predloocv),
        statModel = as.data.frame(stats_model),
        statLOOCV = as.data.frame(stats_loocv)
      )
    } else {
      results[[yy]] <- (list(method = method, model = model))
    }
  }

  if (length(dim(y)) == 1) {
    return(results[[1]])
  }

  return(results)
}
