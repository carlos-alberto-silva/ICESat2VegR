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



#' Fit a Random Forest with optional resampling, tuning, and progress bars
#'
#' `fit_model()` trains a Random Forest regression model and optionally
#' evaluates it via LOOCV, K-fold CV, a train/test split, or bootstrap OOB
#' estimation. It can also tune core RF hyperparameters and shows progress
#' bars for long-running loops.
#'
#' @section Workflow:
#' 1. Optionally tune RF hyperparameters on the *current training subset*
#'    (or on full data when `test$method = "none"`).
#' 2. Always fit a model on the full dataset (`model` and `fitted_full`).
#' 3. If a testing method is requested, refit on resampled training folds and
#'    compute out-of-fold predictions to summarize generalization performance.
#'
#' @param x A data.frame of predictors (rows = samples, cols = features).
#' @param y A numeric vector of responses, `length(y) == nrow(x)`.
#'
#' @param rf_args A named list of base Random Forest arguments used for all
#'   fits (full-data fit and each resample). Elements:
#'   \itemize{
#'     \item `ntree` (integer, default `500`): number of trees.
#'     \item `mtry` (integer or `NULL`): number of variables sampled at each split.
#'           If `NULL`, a safe default `floor(sqrt(p))` is used.
#'     \item `nodesize` (integer, default `5`): minimum terminal node size.
#'     \item `sampsize` (integer, fraction in `(0,1]`, or `NULL`):
#'           sample size per tree. If a single fraction, it is multiplied by the
#'           current training size and rounded. If `NULL`, the package default is used.
#'   }
#'
#' @param test A named list configuring evaluation:
#'   \itemize{
#'     \item `method` (character): one of `none`, `loocv`, `k-fold` (also accepts `kfold`/`k fold`),
#'           `split`, or `bootstrap`.
#'     \item `k` (integer, default `5`): number of folds for K-fold CV.
#'     \item `folds` (list or `NULL`): custom index list of test-fold indices.
#'           If supplied, overrides `k`.
#'     \item `test_size` (numeric in `(0,1)`, default `0.3`): test fraction for
#'           train/test split.
#'     \item `iterations` (integer, default `200`): number of bootstrap resamples.
#'     \item `correction` (logical, default `FALSE`): if `TRUE`, apply the
#'           .632 correction to the RMSE for the bootstrap estimate.
#'     \item `seed` (integer or `NULL`): RNG seed for reproducibility (folds,
#'           split, random tuning).
#'   }
#'
#' @param tune A named list configuring hyperparameter search on the *training*
#'   subset for each fit:
#'   \itemize{
#'     \item `enable` (logical, default `FALSE`): turn tuning on/off.
#'     \item `search` (character, default `grid`): `grid` or `random`.
#'     \item `grid` (data.frame or `NULL`): explicit grid with columns
#'           `mtry`, `ntree`, `nodesize`, `sampsize`. If `NULL`, a sensible
#'           default grid is generated from `p = ncol(x)`.
#'     \item `n_random` (integer, default `20`): number of random draws from
#'           the grid when `search = "random"`.
#'     \item `seed` (integer or `NULL`): RNG seed for the tuning subset draw order.
#'   }
#'
#' @param verbose Logical (default `TRUE`): show progress bars/messages for
#'   tuning, LOOCV, K-fold, and bootstrap loops. Set `FALSE` to silence.
#'
#' @param list_test_models Logical (default `TRUE`): if `TRUE`, save the
#'   *per-resample* fitted model objects as a list in `models_test`.
#'   This can be large, especially for `bootstrap` with many iterations.
#'
#' @details
#' Tuning minimizes OOB MSE for each candidate configuration on the current
#' training subset and then refits the model using the best settings.
#' When `test$method = "bootstrap"` and `correction = TRUE`, the .632 corrected
#' RMSE is reported while keeping other OOB statistics unchanged.
#'
#' @return A list with:
#' \describe{
#'   \item{method}{`rf`}
#'   \item{rf_args}{Final RF arguments used for the full-data fit (after tuning).}
#'   \item{tune_table}{(data.frame or `NULL`) tuning results sorted by OOB MSE.}
#'   \item{test}{Echo of `test` configuration (with resolved values).}
#'   \item{model}{`randomForest` object fit on the full dataset (or train subset for `split`).}
#'   \item{fitted_full}{Numeric vector of in-sample predictions from the full-data fit.}
#'   \item{stats_train}{data.frame of training statistics (RMSE, Bias, %RMSE, %Bias, r, r2).}
#'   \item{stats_test}{(data.frame or `NULL`) test-set or OOF statistics depending on method.}
#'   \item{loocv_pred}{(numeric) LOOCV predictions (for `method = "loocv"`).}
#'   \item{cv_pred}{(numeric) K-fold OOF predictions (for `method = "k-fold"`).}
#'   \item{train_index,test_index}{(integer) indices for `method = "split"`.}
#'   \item{pred_train,pred_test}{(numeric) predictions for train/test in `split`.}
#'   \item{oob_pred}{(numeric) OOB mean prediction per observation in `bootstrap`.}
#'   \item{models_test}{(list or `NULL`) models fitted per resample/fold/iteration when `list_test_models=TRUE`.}
#' }
#'
#' @section Progress Bars:
#' Uses `utils::txtProgressBar()`; disable with `verbose = FALSE`. The internal
#' helper is lightweight and has no external dependencies.
#'
#' @examples
#' set.seed(42)
#' n <- 200
#' x <- data.frame(NDVI = runif(n, 0.2, 0.9),
#'                 EVI  = runif(n, 0.1, 0.8),
#'                 NBR  = runif(n, -0.5, 0.9),
#'                 SLP  = runif(n, 0, 30))
#' y <- with(x, 5 + 20*NDVI + 10*EVI^1.5 - 0.05*SLP + rnorm(n, 0, 2))
#'
#' # Full-data fit (no resampling)
#' fit_none <- fit_model(x, y, rf_args = list(ntree = 400, mtry = 2))
#' fit_none$stats_train
#'
#' # LOOCV
#' \dontshow{
#'   # Make it smaller for CRAN checks
#'   x <- x[1:20, ]
#'   y <- y[1:20]
#' }
#' fit_loocv <- fit_model(x, y, test = list(method = "loocv"))
#' fit_loocv$stats_test
#'
#' # 5-fold CV
#' fit_k_fold <- fit_model(x, y, test = list(method = "k-fold", k = 5, seed = 42))
#' fit_k_fold$stats_test
#'
#' # Train/Test split
#' fit_split <- fit_model(x, y, test = list(method = "split", test_size = 0.25, seed = 42))
#' fit_split$stats_test
#'
#' # Bootstrap (with tuning and .632 correction)
#' ntree <- 400
#' iterations <- 300
#' \dontshow{
#'   # Smaller tree and iterations for checks
#'   ntree <- 50
#'   iterations <- 10
#' }
#' fit_boot <- fit_model(
#'   x, y,
#'   rf_args = list(ntree = ntree),
#'   test    = list(method = "bootstrap", iterations = iterations, correction = TRUE),
#'   tune    = list(enable = TRUE, search = "random", n_random = 12)
#' )
#'
#' # K-fold CV with saved models
#' fit_k <- fit_model(x, y, test = list(method = "k-fold", k = 5, seed = 42), list_test_models = TRUE)
#' length(fit_k$models_test)  # one model per fold
#'
#' @seealso [randomForest::randomForest()], [stats::lm()], [stats::cor()]
#' @importFrom stats lm cor
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom graphics plot
#' @export
fit_model <- function(
    x, y,
    rf_args = list(ntree = 500, mtry = NULL, nodesize = 5, sampsize = NULL),
    test    = list(method = "none", k = 5, test_size = 0.3, seed = NULL, folds = NULL,
                   iterations = 200, correction = FALSE),
    tune    = list(enable = FALSE, search = "grid", grid = NULL, n_random = 20, seed = NULL),
    verbose = TRUE,
    list_test_models = TRUE
) {
  stopifnot(is.data.frame(x), is.numeric(y), nrow(x) == length(y))
  if (!requireNamespace("randomForest", quietly = TRUE))
    stop("Package 'randomForest' is required.")
  n <- nrow(x); p <- ncol(x)

  # ---------- helper ----------
  `%||%` <- function(a, b) if (!is.null(a)) a else b

  .new_pb <- function(n_steps, enabled = TRUE, label = NULL) {
    if (!enabled || n_steps <= 1L) {
      return(list(tick = function() NULL, close = function() NULL))
    }
    if (!is.null(label)) cat(sprintf("%s (%d steps)\n", label, n_steps))
    pb <- utils::txtProgressBar(min = 0, max = n_steps, style = 3)
    i  <- 0L
    list(
      tick  = function() { i <<- i + 1L; utils::setTxtProgressBar(pb, i) },
      close = function() { close(pb); cat("\n") }
    )
  }

  # normalize RF args
  norm_rf_args <- function(args, n_train) {
    a <- args
    # default mtry if missing or NA
    if (is.null(a$mtry) || is.na(a$mtry)) {
      a$mtry <- max(1L, min(floor(sqrt(p)), p))
    }
    if (!is.null(a$mtry)) a$mtry <- max(1L, min(as.integer(a$mtry), p))
    if (!is.null(a$nodesize)) a$nodesize <- as.integer(a$nodesize)

    # sampsize handling:
    # - NULL or NA => let randomForest use its default
    # - fraction in (0,1] => scale by n_train
    # - integer => clamp to [1, n_train]
    if (is.null(a$sampsize) || (length(a$sampsize) == 1L && is.na(a$sampsize))) {
      a$sampsize <- NULL
    } else {
      if (length(a$sampsize) == 1L && is.numeric(a$sampsize) &&
          is.finite(a$sampsize) && a$sampsize > 0 && a$sampsize <= 1) {
        a$sampsize <- as.integer(round(a$sampsize * n_train))
      } else {
        a$sampsize <- as.integer(a$sampsize)
      }
      if (!length(a$sampsize) || any(!is.finite(a$sampsize)) || any(is.na(a$sampsize))) {
        a$sampsize <- NULL
      } else {
        a$sampsize <- pmin(pmax(a$sampsize, 1L), n_train)
      }
    }
    a
  }

  # tuning helpers (train-only)
  default_grid <- function(n_train) {
    mtry_vals     <- unique(pmax(1L, pmin(p, c(ceiling(sqrt(p)), ceiling(p/3), p))))
    ntree_vals    <- c(200L, 400L, 600L)
    nodesize_vals <- c(1L, 5L, 10L)
    sampsize_vals <- c(0.632, 0.8, 1.0)
    expand.grid(mtry = mtry_vals, ntree = ntree_vals,
                nodesize = nodesize_vals, sampsize = sampsize_vals,
                KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  }
  oob_mse <- function(fit) if (!is.null(fit$mse)) utils::tail(fit$mse, 1) else NA_real_

  tune_rf <- function(xtr, ytr, tune_cfg) {
    if (!isTRUE(tune_cfg$enable)) return(list(best = list(), table = NULL))
    n_tr <- nrow(xtr)
    grid <- tune_cfg$grid; if (is.null(grid)) grid <- default_grid(n_tr)
    if (identical(tune_cfg$search, "random")) {
      set.seed(tune_cfg$seed %||% NULL)
      if (nrow(grid) > (tune_cfg$n_random %||% 20))
        grid <- grid[sample.int(nrow(grid), tune_cfg$n_random), , drop = FALSE]
    }
    pb <- .new_pb(nrow(grid), enabled = verbose, label = "Tuning randomForest")
    on.exit(pb$close(), add = TRUE)
    res <- vector("list", nrow(grid))
    for (i in seq_len(nrow(grid))) {
      g <- as.list(grid[i, , drop = FALSE])
      if (!is.null(g$sampsize) && g$sampsize > 0 && g$sampsize < 1)
        g$sampsize <- as.integer(round(g$sampsize * n_tr))
      fit <- do.call(randomForest::randomForest, c(list(x = xtr, y = ytr), g))
      res[[i]] <- data.frame(mtry=g$mtry, ntree=g$ntree, nodesize=g$nodesize,
                             sampsize=g$sampsize, oob_mse=oob_mse(fit), stringsAsFactors=FALSE)
      pb$tick()
    }
    tab <- do.call(rbind, res)
    tab <- tab[order(tab$oob_mse), , drop = FALSE]
    best <- as.list(tab[1, c("mtry","ntree","nodesize","sampsize")])
    list(best = best, table = tab)
  }

  fit_train <- function(xtr, ytr, base_args, tune_cfg) {
    tun <- tune_rf(xtr, ytr, tune_cfg)
    a   <- base_args
    if (length(tun$best)) {
      a$mtry     <- tun$best$mtry
      a$ntree    <- tun$best$ntree
      a$nodesize <- tun$best$nodesize
      a$sampsize <- tun$best$sampsize
    }
    a <- norm_rf_args(a, n_train = nrow(xtr))
    fit <- do.call(randomForest::randomForest, c(list(x = xtr, y = ytr), a))
    list(fit = fit, args = a, tune_table = tun$table)
  }

  # ---- method normalization & seed ----
  method_raw <- tolower(test$method %||% "none")
  method_key <- gsub("[^a-z]", "", method_raw)  # "k-fold"/"k fold" -> "kfold"
  if (!is.null(test$seed)) set.seed(test$seed)

  # ---- full-data fit (always) ----
  full_pack    <- fit_train(x, y, rf_args, utils::modifyList(tune, list(enable = tune$enable)))
  full_fit     <- full_pack$fit
  full_args    <- full_pack$args
  fitted_full  <- as.numeric(stats::predict(full_fit, newdata = x))
  stats_train  <- fit_metrics(y, fitted_full)

  # no test
  if (method_key == "none") {
    return(list(
      method      = "rf",
      rf_args     = full_args,
      tune_table  = full_pack$tune_table,
      test        = list(method = "none"),
      model       = full_fit,
      fitted_full = fitted_full,
      stats_train = stats_train,
      stats_test  = NULL,
      models_test = NULL
    ))
  }

  # ---- LOOCV ----
  if (method_key == "loocv") {
    pred_oof <- numeric(n)
    models_list <- if (list_test_models) vector("list", n) else NULL
    pb <- .new_pb(n, enabled = verbose, label = "LOOCV")
    on.exit(pb$close(), add = TRUE)
    for (i in seq_len(n)) {
      idx_tr <- setdiff(seq_len(n), i)
      pack   <- fit_train(x[idx_tr, , drop = FALSE], y[idx_tr], rf_args, tune)
      pred_oof[i] <- as.numeric(stats::predict(pack$fit, newdata = x[i, , drop = FALSE]))
      if (list_test_models) models_list[[i]] <- pack$fit
      pb$tick()
    }
    return(list(
      method      = "rf",
      rf_args     = full_args,
      tune_table  = full_pack$tune_table,
      test        = list(method = "loocv"),
      model       = full_fit,
      fitted_full = fitted_full,
      loocv_pred  = pred_oof,
      stats_train = stats_train,
      stats_test  = fit_metrics(y, pred_oof),
      models_test = models_list
    ))
  }

  # ---- K-fold ----
  if (method_key == "kfold") {
    k <- as.integer(test$k %||% 5L)
    if (!is.null(test$folds)) {
      folds <- test$folds
      if (!is.list(folds)) stop("'folds' must be a list of index vectors.")
    } else {
      ord <- sample.int(n)
      folds <- split(ord, rep(1:k, length.out = n))
    }
    pred_oof <- rep(NA_real_, n)
    models_list <- if (list_test_models) vector("list", length(folds)) else NULL
    pb <- .new_pb(length(folds), enabled = verbose, label = sprintf("K-fold CV (k=%d)", length(folds)))
    on.exit(pb$close(), add = TRUE)
    for (f in seq_along(folds)) {
      idx_te <- folds[[f]]
      idx_tr <- setdiff(seq_len(n), idx_te)
      pack   <- fit_train(x[idx_tr, , drop = FALSE], y[idx_tr], rf_args, tune)
      pred_oof[idx_te] <- as.numeric(stats::predict(pack$fit, newdata = x[idx_te, , drop = FALSE]))
      if (list_test_models) models_list[[f]] <- pack$fit
      pb$tick()
    }
    return(list(
      method      = "rf",
      rf_args     = full_args,
      tune_table  = full_pack$tune_table,
      test        = list(method = "k-fold", k = length(folds), folds = folds),
      model       = full_fit,
      fitted_full = fitted_full,
      cv_pred     = pred_oof,
      stats_train = stats_train,
      stats_test  = fit_metrics(y, pred_oof),
      models_test = models_list
    ))
  }

  # ---- Split ----
  if (method_key == "split") {
    test_size <- as.numeric(test$test_size %||% 0.3)
    test_size <- max(0.05, min(0.95, test_size))
    n_te   <- max(1L, as.integer(round(n * test_size)))
    idx_te <- sample(seq_len(n), size = n_te)
    idx_tr <- setdiff(seq_len(n), idx_te)

    pack    <- fit_train(x[idx_tr, , drop = FALSE], y[idx_tr], rf_args, tune)
    fit     <- pack$fit
    pred_tr <- as.numeric(stats::predict(fit, newdata = x[idx_tr, , drop = FALSE]))
    pred_te <- as.numeric(stats::predict(fit, newdata = x[idx_te, , drop = FALSE]))
    models_list <- if (list_test_models) list(train_model = fit) else NULL

    return(list(
      method       = "rf",
      rf_args      = pack$args,
      tune_table   = pack$tune_table,
      test         = list(method = "split", test_size = test_size, seed = test$seed),
      model        = fit,             # trained on train subset
      fitted_full  = fitted_full,     # full-data reference
      train_index  = idx_tr,
      test_index   = idx_te,
      pred_train   = pred_tr,
      pred_test    = pred_te,
      stats_train  = fit_metrics(y[idx_tr], pred_tr),
      stats_test   = fit_metrics(y[idx_te], pred_te),
      models_test  = models_list
    ))
  }

  # ---- Bootstrap ----
  if (method_key == "bootstrap") {
    iterations     <- as.integer(test$iterations %||% 200L)
    use_correction <- isTRUE(test$correction)

    oof_pred <- vector("list", n)  # per-row OOB predictions across resamples
    models_list <- if (list_test_models) vector("list", iterations) else NULL
    pb <- .new_pb(iterations, enabled = verbose, label = sprintf("Bootstrap (%d iterations)", iterations))
    on.exit(pb$close(), add = TRUE)
    for (b in seq_len(iterations)) {
      idx_tr <- sample.int(n, replace = TRUE)
      inbag  <- logical(n); inbag[idx_tr] <- TRUE
      idx_oob <- which(!inbag)
      if (length(idx_oob) > 0) {
        pack <- fit_train(x[idx_tr, , drop = FALSE], y[idx_tr], rf_args, tune)
        preds <- as.numeric(stats::predict(pack$fit, newdata = x[idx_oob, , drop = FALSE]))
        for (j in seq_along(idx_oob)) {
          i <- idx_oob[j]
          oof_pred[[i]] <- c(oof_pred[[i]], preds[j])
        }
        if (list_test_models) models_list[[b]] <- pack$fit
      } else if (list_test_models) {
        models_list[[b]] <- NULL
      }
      pb$tick()
    }
    pred_oob <- rep(NA_real_, n)
    for (i in seq_len(n)) if (length(oof_pred[[i]])) pred_oob[i] <- mean(oof_pred[[i]])

    if (!use_correction) {
      keep <- which(is.finite(pred_oob))
      stats_test <- fit_metrics(y[keep], pred_oob[keep])
    } else {
      keep   <- which(is.finite(pred_oob))
      s_tr   <- fit_metrics(y, fitted_full)
      s_oob  <- fit_metrics(y[keep], pred_oob[keep])
      rmse_632 <- 0.368 * s_tr$value[s_tr$stat=="rmse"] +
        0.632 * s_oob$value[s_oob$stat=="rmse"]
      stats_test <- s_oob
      stats_test$value[stats_test$stat=="rmse"] <- as.numeric(rmse_632)
      stats_test$note <- if ("note" %in% names(stats_test)) stats_test$note else NA
      stats_test$note[stats_test$stat=="rmse"] <- "0.632-corrected"
    }

    return(list(
      method       = "rf",
      rf_args      = full_args,
      tune_table   = full_pack$tune_table,
      test         = list(method = "bootstrap", iterations = iterations, correction = use_correction),
      model        = full_fit,      # full-data model for inference
      fitted_full  = fitted_full,
      oob_pred     = pred_oob,      # NA if never OOB
      stats_train  = stats_train,
      stats_test   = stats_test,
      models_test  = models_list
    ))
  }

  stop("Unknown test$method: use one of 'none','loocv','k-fold','split','bootstrap' (hyphen/space accepted).")
}

