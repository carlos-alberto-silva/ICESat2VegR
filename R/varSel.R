## ============================================================================
##  Internal importance scaling helper (Breiman / randomForest only)
## ============================================================================

#' @title Internal scaling of Random Forest importance values
#' @description
#' Internal helper for scaling permutation-based **randomForest** importance
#' values. Supports the `mir` and `se`` scaling options.
#'
#' @param x A [randomForest::randomForest] object (classification
#'   or regression) with permutation importance computed.
#' @param scaling Character; type of importance scaling, either `mir` or
#'   `se`.
#' @param sort Logical; if `TRUE`, the output is sorted by importance.
#'
#' @return
#' A data.frame with:
#'
#' - `parameter` - variable name
#' - `importance` - scaled importance value
#'
#' @details
#' The scaled importance measures are calculated as:
#' \deqn{\mathrm{mir}_j = \frac{i_j}{\max(i)}}
#' \deqn{\mathrm{se}_j = \frac{i_j / \mathrm{SE}_j}{\sum_k (i_k / \mathrm{SE}_k)}}
#'
#' @keywords internal
.rfImpScale <- function(x,
                        scaling = c("mir", "se"),
                        sort    = FALSE) {

  scaling <- match.arg(scaling)

  if (!inherits(x, "randomForest")) {
    stop("x must be a 'randomForest' object.")
  }

  mtype <- tolower(x$type)

  ## ---- Extract raw importance ----
  if (mtype == "regression") {
    if (!"%IncMSE" %in% colnames(x$importance)) {
      stop("randomForest object does not contain '%IncMSE' importance.")
    }
    rf.imp <- x$importance[, "%IncMSE"]

  } else if (mtype %in% c("classification", "unsupervised", "probability estimation")) {
    if (!"MeanDecreaseAccuracy" %in% colnames(x$importance)) {
      stop("randomForest object does not contain 'MeanDecreaseAccuracy' importance.")
    }
    rf.imp <- x$importance[, "MeanDecreaseAccuracy"]

  } else {
    stop("Unsupported randomForest type: ", mtype)
  }

  ## ---- Scaling ----
  if (scaling == "mir") {

    i <- rf.imp / max(rf.imp, na.rm = TRUE)

  } else if (scaling == "se") {

    if (mtype == "regression") {
      rf.impSD <- x$importanceSD
    } else if (mtype == "classification") {
      rf.impSD <- x$importanceSD[, "MeanDecreaseAccuracy"]
    } else {
      stop("Unsupported randomForest type for 'se' scaling: ", mtype)
    }

    rf.impSD[rf.impSD == 0] <- 1e-9
    i <- (rf.imp / rf.impSD) / sum(rf.imp / rf.impSD, na.rm = TRUE)
  }

  out <- data.frame(
    parameter  = names(i),
    importance = as.numeric(i),
    row.names  = NULL
  )

  if (sort) {
    out <- out[order(out$importance), ]
    rownames(out) <- NULL
  }

  out
}


## ============================================================================
##  Internal confusion-matrix accuracy helper (no rfUtilities dependency)
## ============================================================================

#' @title Internal confusion matrix accuracy
#' @description
#' Internal helper to compute overall PCC, producer's accuracy and kappa from
#' a confusion matrix. Adapted conceptually from `rfUtilities::accuracy`.
#'
#' @param cmat Confusion matrix with rows = reference (true) classes and
#'   columns = predicted classes.
#'
#' @return
#' A list with:
#' - `PCC` - overall percent correctly classified.
#' - `kappa` - Cohen's kappa coefficient.
#' - `producers.accuracy` - vector of producer's accuracy per class.
#'
#'
#' @keywords internal
.rfAccuracy <- function(cmat) {

  if (!is.matrix(cmat) && !is.table(cmat)) {
    stop("'cmat' must be a matrix or table.")
  }

  cmat <- as.matrix(cmat)

  N         <- sum(cmat)
  correct   <- sum(diag(cmat))
  row.tot   <- rowSums(cmat)
  col.tot   <- colSums(cmat)

  # Overall accuracy (PCC)
  PCC <- (correct / N) * 100

  # Producer's accuracy (per class)
  producers <- diag(cmat) / row.tot * 100
  names(producers) <- rownames(cmat)

  # Kappa
  Po <- correct / N
  Pe <- sum(row.tot * col.tot) / (N^2)
  kappa <- (Po - Pe) / (1 - Pe)

  list(
    PCC                = PCC,
    kappa              = kappa,
    producers.accuracy = producers
  )
}


## ============================================================================
##  varSel: Random Forest variable selection (Breiman-only)
## ============================================================================

#' @title Random Forest Variable Selection (Breiman-only)
#' @description
#' Implements the Random Forest model selection approach of Murphy et al. (2010),
#' using Breiman's original **randomForest** implementation. This is a simplified
#' adaptation of `rfUtilities::rf.modelSel`, restricted to the Breiman
#' implementation and modified for the ICESat2VegR package.
#'
#' It returns the selected variables and importance metrics, but does not fit
#' a final model.
#'
#' @param xdata Matrix or data.frame of predictor variables (columns = predictors).
#' @param ydata Response vector. For classification, `ydata` must be a factor;
#'   otherwise the model is fit in regression mode.
#' @param imp.scale Character; type of scaling for importance values, either
#'   `mir` or `se`. Default is `mir`.
#' @param r Numeric vector of importance percentiles to test, e.g.,
#'   `c(0.25, 0.50, 0.75)`. These percentiles are used to define thresholds
#'   for building nested models during selection.
#' @param min.imp Optional numeric in \eqn{[ 0,1 ]} used as a minimum scaled
#'   importance cutoff (in MIR or SE scale) to filter the final set of variables.
#'   Variables whose scaled importance is below this threshold are dropped from
#'   the selected set `selvars`. If `NULL` (default), no cutoff is
#'   applied after the Murphy-style model selection.
#' @param seed Optional integer; sets the random seed in the global R environment.
#'   This is strongly recommended for reproducibility.
#' @param parsimony Numeric in (0,1); threshold for selecting among competing
#'   models. If specified, models whose errors are within `parsimony` of
#'   the best model (OOB / class error for classification, variance explained /
#'   MSE for regression) are considered, and the one with the fewest parameters
#'   is chosen.
#' @param kappa Logical; use the chance-corrected \eqn{\kappa} statistic as the
#'   primary optimization criterion for classification instead of percent
#'   correctly classified (PCC). This corrects PCC for random agreement.
#' @param ... Additional arguments passed to [randomForest::randomForest()],
#'   e.g. `ntree = 1000`, `replace = TRUE`, `proximity = TRUE`.
#'
#' @return An object of class `varSel` (a list) with components:
#'   - `selvars` Character vector of selected variable names (after
#'     applying `min.imp`, if provided).
#'   - `test` Data.frame of model-selection diagnostics, containing
#'     error metrics, threshold, number of parameters, and the variables used
#'     in each candidate model (for inspection only).
#'   - `importance` Data.frame of scaled importance values for all
#'     variables in the full model (columns `parameter`, `importance`).
#'   - `sel.importance` Data.frame of scaled importance values for the
#'     selected variables.
#'   - `parameters` List of variables used in each candidate model.
#'   - `scaling` Character; the importance scaling used (`mir`
#'     or `se`).
#'
#'
#' @details
#' For classification, ensure that `ydata` is a factor; otherwise the model
#' is fit in regression mode.
#'
#' **Selection strategy:**
#'
#'   - For classification, candidate models are compared using OOB PCC (or
#'     \eqn{\kappa}, if `kappa = TRUE`) and maximum within-class error,
#'     with preference for more parsimonious models.
#'   - For regression, candidate models are compared using percent variance
#'     explained and MSE, with preference for more parsimonious models.
#'   - After the best model is chosen, the optional `min.imp` cutoff is
#'     applied to the scaled importance (MIR/SE) from the full model to remove
#'     weak predictors from the final set `selvars`.
#'
#'
#' Typical choices for `min.imp`:
#'   - MIR: `min.imp` in \eqn{[0.1, 0.3]} (e.g., keep variables with
#'     at least 20\% of the maximum importance).
#'   - SE: `min.imp` in \eqn{[0.01, 0.05]}, remembering that SE-scaled
#'     importance values sum to 1.
#'
#'
#' @examples
#' require(randomForest)
#'
#' data(airquality)
#' airquality <- na.omit(airquality)
#'
#' xdata <- airquality[, 2:6]
#' ydata <- airquality[, 1]
#' ntree <- 500
#' \dontshow{
#'   # Lightweight version for R CMD check
#'   xdata <- airquality[1:50, 2:6]
#'   ydata <- airquality[1:50, 1]
#'   ntree <- 50
#' }
#' ## Regression example with MIR scaling and importance cutoff
#' vs.regress <- varSel(
#'   xdata     = xdata,
#'   ydata     = ydata,
#'   imp.scale = "mir",
#'   ntree     = ntree,
#'   min.imp   = 0.2
#' )
#'
#' vs.regress$selvars
#'
#' ## Plot all variables, highlighting selected ones
#' plot(vs.regress, which = "importance")
#'
#' @seealso
#' [randomForest::randomForest()] and [plot.varSel()].
#'
#' @export
varSel <- function(xdata,
                   ydata,
                   imp.scale = c("mir", "se"),
                   r         = c(0.25, 0.50, 0.75),
                   min.imp   = NULL,
                   seed      = NULL,
                   parsimony = NULL,
                   kappa     = FALSE,
                   ...) {

  if (missing(ydata) || missing(xdata)) {
    stop("Both 'ydata' and 'xdata' must be provided.")
  }

  imp.scale <- match.arg(imp.scale)

  if (!is.null(seed)) set.seed(seed)

  # Percentile thresholds
  r         <- unique(c(0, r, 1))
  expected  <- seq(0, 1, 0.01)
  range.idx <- findInterval(expected, r, all.inside = TRUE)

  # Build argument list for randomForest
  dots <- as.list(match.call(expand.dots = TRUE))[-1]
  rm.idx <- names(dots) %in% c(
    "xdata", "ydata", "imp.scale", "r",
    "min.imp", "seed", "parsimony", "kappa"
  )
  if (any(rm.idx)) dots <- dots[!rm.idx]

  if (!"x" %in% names(dots)) dots[["x"]] <- xdata
  if (!"y" %in% names(dots)) dots[["y"]] <- ydata
  if (!"importance" %in% names(dots)) dots[["importance"]] <- TRUE

  RFtype <- is.factor(ydata)  # TRUE = classification, FALSE = regression

  model.vars      <- list()
  imp             <- NULL
  sel.importance  <- NULL
  selvars         <- NULL

  ## ---------------------------------------------------------------------------
  ## CLASSIFICATION
  ## ---------------------------------------------------------------------------
  if (RFtype) {

    ln <- 0

    rf.all <- do.call(randomForest::randomForest, dots)
    model.vars[[ln <- ln + 1]] <- rownames(rf.all$importance)

    cm <- .rfAccuracy(rf.all$confusion[, 1:(ncol(rf.all$confusion) - 1)])

    if (kappa) {
      class.errors <- data.frame(
        1 - t(c(cm$kappa, (cm$producers.accuracy / 100)))
      )
      names(class.errors)[1] <- "kappa"
    } else {
      class.errors <- data.frame(
        1 - t(c(cm$PCC, cm$producers.accuracy / 100))
      )
      names(class.errors)[1] <- "pcc"
    }

    errors <- data.frame(
      error       = class.errors[, 1],
      class.error = max(class.errors[, 2:ncol(class.errors)]),
      threshold   = 1,
      nparameters = ncol(xdata)
    )

    nan.idx <- which(is.nan(errors))
    if (length(nan.idx) > 0) errors[1, ][nan.idx] <- 1
    inf.idx <- which(is.infinite(errors))
    if (length(inf.idx) > 0) errors[1, ][inf.idx] <- 0

    # Importance scaling (full model)
    imp <- .rfImpScale(rf.all, scaling = imp.scale, sort = FALSE)
    imp <- imp[order(imp$importance), ]
    imp$p <- findInterval(imp$importance, r, all.inside = TRUE)

    # Warnings about ranges without variables
    miss.idx <- which(
      !sort(unique(findInterval(expected, r, all.inside = TRUE))) %in%
        sort(unique(imp$p))
    )
    if (length(miss.idx) > 0) {
      miss.range <- tapply(
        expected,
        findInterval(expected, r, all.inside = TRUE),
        range
      )[miss.idx]

      miss.range <- paste0(
        "Missing importance ranges: ",
        unlist(lapply(miss.range, function(x) paste0(x[1], "-", x[2])))
      )
      for (i in miss.range) message(i)
    }

    # Sequencial models by importance threshold
    for (p in sort(unique(imp$p))) {
      idx.p     <- which(imp$p == p)
      sel.imp.p <- imp[idx.p, , drop = FALSE]

      if (!exists("sel.vars")) {
        sel.vars <- imp$parameter
      }

      # Skip if parameters don't change
      if (!any((sel.vars %in% sel.imp.p$parameter) == FALSE)) {
        rng <- range(expected[range.idx %in% p])
        warning(sprintf(
          "The %.2f-%.2f threshold has <= 1 parameter and cannot be evaluated",
          rng[1], rng[2]
        ))
        next
      }

      sel.vars <- sel.imp.p$parameter
      np <- length(sel.vars)

      if (np > 1) {
        thres <- p

        dots[["x"]] <- xdata[, sel.vars, drop = FALSE]
        rf.model    <- do.call(randomForest::randomForest, dots)
        model.vars[[ln <- ln + 1]] <- rownames(rf.model$importance)

        cm <- .rfAccuracy(rf.model$confusion[, 1:(ncol(rf.model$confusion) - 1)])

        if (kappa) {
          e <- data.frame(
            1 - t(c(cm$kappa, (cm$producers.accuracy / 100)))
          )
          names(e)[1] <- "kappa"
        } else {
          e <- data.frame(
            1 - t(c(cm$PCC, cm$producers.accuracy / 100))
          )
          names(e) <- c("pcc", names(cm$producers.accuracy))
        }

        e <- data.frame(
          error       = e[, 1],
          class.error = max(e[, 2:ncol(e)]),
          threshold   = thres,
          nparameters = np
        )

        nan.idx <- which(is.nan(e))
        if (length(nan.idx) > 0) e[1, ][nan.idx] <- 1
        inf.idx <- which(is.infinite(e))
        if (length(inf.idx) > 0) e[1, ][inf.idx] <- 0

        errors <- rbind(errors, e)
      } else {
        rng <- range(expected[range.idx %in% p])
        warning(sprintf(
          "The %.2f-%.2f threshold has <= 1 parameter and cannot be evaluated",
          rng[1], rng[2]
        ))
      }
    }

    # Stacks model.vars and data.frame for inspecting
    n <- max(unlist(lapply(model.vars, length)))
    for (l in seq_along(model.vars)) {
      x <- model.vars[[l]]
      length(x) <- n
      model.vars[[l]] <- x
    }
    model.vars.df <- as.data.frame(do.call(rbind, model.vars), stringsAsFactors = FALSE)
    names(model.vars.df) <- paste0("parameter", 1:ncol(model.vars.df))

    errors <- data.frame(errors, model.vars.df)
    rownames(errors) <- seq_len(nrow(errors))
    errors <- errors[order(errors$class.error, errors$error, errors$nparameters), ]

    # Parsimony threshold
    if (!is.null(parsimony)) {
      if (parsimony < 1e-8 || parsimony > 0.9) {
        stop("'parsimony' must range between 0 and 1.")
      }
      oob <- "TRUE"
      for (i in 2:nrow(errors)) {
        if (abs((errors[i, 1] - errors[1, 1]) / errors[1, 2]) <= parsimony &&
            abs((errors[i, 2] - errors[1, 2]) / errors[1, 3]) <= parsimony) {
          oob <- append(oob, "TRUE")
        } else {
          oob <- append(oob, "FALSE")
        }
      }
      final.df <- errors[oob == "TRUE", ]
      final.th <- min(final.df[final.df$nparameters == min(final.df$nparameters), ]$threshold)
    } else {
      final.th <- errors$threshold[1]
    }

    # Variables in the final model
    sel.row    <- which(errors$threshold == final.th)[1]
    param.cols <- grep("^parameter", names(errors))
    selvec     <- as.character(unlist(errors[sel.row, param.cols]))
    selvec     <- selvec[!is.na(selvec) & selvec != ""]
    selvars    <- unique(selvec)

    # Put chosen row at the top of 'test'
    errors <- rbind(errors[sel.row, ], errors[-sel.row, ])

    # Importance before scalling (min.imp)
    sel.importance <- imp[imp$parameter %in% selvars, c("parameter", "importance")]

    ## ---------------------------------------------------------------------------
    ## REGRESSION
    ## ---------------------------------------------------------------------------
  } else {

    ln <- 0

    rf.all <- randomForest::randomForest(
      x = xdata,
      y = ydata,
      importance = TRUE,
      ...
    )
    model.vars[[ln <- ln + 1]] <- rownames(rf.all$importance)
    e <- c(
      varexp = stats::median(rf.all$rsq),
      mse    = mean(rf.all$mse)
    )

    errors <- data.frame(
      varexp     = e[1],
      mse        = e[2],
      threshold  = 1,
      nparameters = ncol(xdata)
    )

    # Importance scaling (full model)
    imp <- .rfImpScale(rf.all, scaling = imp.scale, sort = FALSE)
    imp <- imp[order(imp$importance), ]
    imp$p <- findInterval(imp$importance, r, all.inside = TRUE)

    miss.idx <- which(
      !sort(unique(findInterval(expected, r, all.inside = TRUE))) %in%
        sort(unique(imp$p))
    )
    if (length(miss.idx) > 0) {
      miss.range <- tapply(
        expected,
        findInterval(expected, r, all.inside = TRUE),
        range
      )[miss.idx]

      miss.range <- paste0(
        "Missing importance ranges: ",
        unlist(lapply(miss.range, function(x) paste0(x[1], "-", x[2])))
      )
      for (i in miss.range) message(i)
    }

    for (p in sort(unique(imp$p))) {
      idx.p     <- which(imp$p == p)
      sel.imp.p <- imp[idx.p, , drop = FALSE]

      if (!exists("sel.vars")) {
        sel.vars <- imp$parameter
      }

      if (!any((sel.vars %in% sel.imp.p$parameter) == FALSE)) {
        rng <- range(expected[range.idx %in% p])
        warning(sprintf(
          "The %.2f-%.2f threshold has <= 1 parameter and cannot be evaluated",
          rng[1], rng[2]
        ))
        next
      }

      sel.vars <- sel.imp.p$parameter
      np <- length(sel.vars)

      if (np > 1) {
        thres <- p

        rf.model <- randomForest::randomForest(
          x = xdata[, sel.vars, drop = FALSE],
          y = ydata,
          importance = TRUE,
          ...
        )
        model.vars[[ln <- ln + 1]] <- rownames(rf.model$importance)

        e <- c(
          varexp = stats::median(rf.model$rsq),
          mse    = mean(rf.model$mse)
        )

        e <- data.frame(
          varexp     = e[1],
          mse        = e[2],
          threshold  = thres,
          nparameters = np
        )
        errors <- rbind(errors, e)
      } else {
        rng <- range(expected[range.idx %in% p])
        warning(sprintf(
          "The %.2f-%.2f threshold has <= 1 parameter and cannot be evaluated",
          rng[1], rng[2]
        ))
      }
    }

    n <- max(unlist(lapply(model.vars, length)))
    for (l in seq_along(model.vars)) {
      x <- model.vars[[l]]
      length(x) <- n
      model.vars[[l]] <- x
    }
    model.vars.df <- as.data.frame(do.call(rbind, model.vars), stringsAsFactors = FALSE)
    names(model.vars.df) <- paste0("parameter", 1:ncol(model.vars.df))

    errors <- data.frame(errors, model.vars.df)
    rownames(errors) <- seq_len(nrow(errors))
    errors <- errors[order(-errors$varexp, errors$mse, errors$nparameters), ]

    if (!is.null(parsimony)) {
      if (parsimony < 1e-8 || parsimony > 0.9) {
        stop("'parsimony' must range between 0 and 1.")
      }
      oob <- "TRUE"
      for (i in 2:nrow(errors)) {
        if (abs((errors[i, 1] - errors[1, 1]) / errors[1, 1]) <= parsimony &&
            abs((errors[i, 2] - errors[1, 2]) / errors[1, 2]) <= parsimony) {
          oob <- append(oob, "TRUE")
        } else {
          oob <- append(oob, "FALSE")
        }
      }
      final.df <- errors[oob == "TRUE", ]
      final.th <- min(final.df[final.df$nparameters == min(final.df$nparameters), ]$threshold)
    } else {
      final.th <- errors$threshold[1]
    }

    sel.row    <- which(errors$threshold == final.th)[1]
    param.cols <- grep("^parameter", names(errors))
    selvec     <- as.character(unlist(errors[sel.row, param.cols]))
    selvec     <- selvec[!is.na(selvec) & selvec != ""]
    selvars    <- unique(selvec)

    errors <- rbind(errors[sel.row, ], errors[-sel.row, ])

    sel.importance <- imp[imp$parameter %in% selvars, c("parameter", "importance")]
  }

  ## -------------------------------------------------------------------------
  ##  Apply cutoff if specified (min.imp)
  ## -------------------------------------------------------------------------
  if (!is.null(min.imp)) {
    if (!is.numeric(min.imp) || length(min.imp) != 1L ||
        min.imp < 0 || min.imp > 1) {
      stop("'min.imp' must be a single numeric between 0 and 1.")
    }

    # Use full model importance for filtering final set
    keep <- imp$importance >= min.imp
    strong.vars <- imp$parameter[keep]

    selvars <- selvars[selvars %in% strong.vars]

    if (length(selvars) == 0) {
      warning("No variables met the 'min.imp' importance cutoff; 'selvars' is empty.")
    }

    sel.importance <- imp[imp$parameter %in% selvars, c("parameter", "importance")]
  }

  mdl.sel <- list(
    selvars        = selvars,
    test           = errors,
    importance     = imp[, c("parameter", "importance")],
    sel.importance = sel.importance,
    parameters     = model.vars,
    scaling        = imp.scale
  )

  class(mdl.sel) <- c("varSel", "list")
  mdl.sel
}


## ============================================================================
##  plot.varSel: horizontal barplot with selection highlighting
## ============================================================================

#' Plot variable importance for `varSel` objects
#'
#' @description
#' Plot scaled variable importance for a `varSel` object using a horizontal
#' barplot. Variables are ordered so that the most important metrics appear at
#' the **top** (largest bars).
#'
#' When `which = "importance"` all variables are shown, and those selected
#' by [varSel()] are highlighted in a different color. An optional
#' color palette (e.g., `viridis::inferno`) can be supplied to color the
#' bars by importance.
#'
#' @param x An object of class `varSel`.
#' @param which Character; either `sel.importance` (default) to plot the
#'   importance of selected variables only, or `importance` to plot the
#'   full-model importance (all predictors).
#' @param main Plot title.
#' @param xlab Label for the x-axis (importance scale).
#' @param col.selected Base color for selected variables when no palette is given.
#' @param col.other Base color for non-selected variables when no palette is given.
#' @param palette Optional color palette. Can be:
#'     - a function `f(n)` returning `n` colors
#'       (e.g., `viridis::inferno`),
#'     - or a character vector of colors used to build a gradient via
#'       [grDevices::colorRampPalette].
#' 
#'   The palette is applied to all bars and then the color of selected variables
#'   is overridden by `col.selected` when `which = "importance"`.
#' @param legend.loc Legend location, e.g. "bottomright"
#' @param ... Additional graphical arguments passed to
#'   [graphics::barplot].
#'
#' @examples
#' require(randomForest)
#'
#' data(airquality)
#' airquality <- na.omit(airquality)
#'
#' xdata <- airquality[, 2:6]
#' ydata <- airquality[, 1]
#'
#' ntree = 200
#' \dontshow{
#'  ntree = 10
#' }
#' vs <- varSel(xdata, ydata, ntree = ntree, min.imp = 0.2)
#'
#' ## Selected variables only
#' plot(vs, which = "sel.importance")
#'
#' ## All variables, highlighting selected ones, with inferno palette
#' if (requireNamespace("viridis", quietly = TRUE)) {
#'   plot(vs, which = "importance", palette = viridis::inferno)
#' }
#'
#' @method plot varSel
#' @importFrom grDevices palette
#' @export
plot.varSel <- function(x,
                        which        = c("sel.importance", "importance"),
                        main         = "Random Forest variable importance",
                        xlab         = "Scaled importance",
                        col.selected = "steelblue",
                        col.other    = "grey80",
                        palette      = NULL,
                        legend.loc = "bottomright",
                        ...) {

  which <- match.arg(which)

  # Base importance table
  if (which == "sel.importance") {
    if (is.null(x$sel.importance)) {
      stop("Object does not contain 'sel.importance'.")
    }
    imp <- x$sel.importance
  } else {
    if (is.null(x$importance)) {
      stop("Object does not contain 'importance'.")
    }
    imp <- x$importance
  }

  if (!all(c("parameter", "importance") %in% names(imp))) {
    stop("Importance data.frame must contain 'parameter' and 'importance' columns.")
  }

  # Sort so that the most important appears at the top:
  # for horizontal bars, the "topmost" bar is the last one drawn,
  # so we sort in increasing order of importance.
  imp <- imp[order(imp$importance, decreasing = FALSE), ]
  rownames(imp) <- NULL

  n <- nrow(imp)

  # Base palette (if provided)
  if (!is.null(palette)) {
    if (is.function(palette)) {
      cols.base <- palette(n)
    } else if (is.character(palette) && length(palette) >= 2) {
      cols.base <- grDevices::colorRampPalette(palette)(n)
    } else {
      stop("'palette' must be a function or a character vector of length >= 2.")
    }
  } else {
    cols.base <- rep(col.other, n)
  }

  # Final colors
  cols <- cols.base

  if (which == "importance") {
    sel <- x$selvars
    if (!is.null(sel) && length(sel) > 0) {
      is.sel <- imp$parameter %in% sel
      cols[is.sel] <- col.selected   # highlight selected ones
    }
  } else {
    # Only selected variables: use col.selected
    cols[] <- col.selected
  }

  old.par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old.par))

  # extra space for labels on the Y axis
  graphics::par(mar = c(5, 10, 4, 2))

  graphics::barplot(
    height    = imp$importance,
    names.arg = imp$parameter,
    horiz     = TRUE,
    las       = 1,     # horizontal labels on y-axis
    col       = cols,
    border    = NA,
    main      = main,
    xlab      = xlab,
    ...
  )

  if (which == "importance") {
    graphics::legend(legend.loc,
           legend = c("Selected", "Not selected"),
           fill   = c(col.selected, col.other),
           bty    = "n")
  }
}
