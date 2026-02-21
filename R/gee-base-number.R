#' Creates an `Earth Engine` server number
#'
#' @param x the number to convert to Earth Engine's number.
#'
#' @return The Earth Engine number
#'
#' @seealso https://developers.google.com/earth-engine/apidocs/ee-number
#' @keywords internal
ee_number <- function(x) {
  UseMethod("ee_number")
}

#' @keywords internal
"ee_number.numeric" <- function(x) if (!is.null(ee)) ee$Number(x)

#' @keywords internal
`+.ee.ee_number.Number` <- function(x, e2) {
  if (!is.null(ee)) {
    if (inherits(x, "numeric")) {
      x <- ee$Number(x)
    }

    x$add(e2)
  }
}

#' @keywords internal
`-.ee.ee_number.Number` <- function(x, e2) {
  if (!is.null(ee)) {
    if (inherits(x, "numeric")) {
      x <- ee$Number(x)
    }

    x$subtract(e2)
  }
}

#' @keywords internal
`*.ee.ee_number.Number` <- function(x, e2) {
  if (!is.null(ee)) {
    if (inherits(x, "numeric")) {
      x <- ee$Number(x)
    }

    x$multiply(e2)
  }
}

#' @keywords internal
`/.ee.ee_number.Number` <- function(x, e2) {
  if (!is.null(ee)) {
    if (inherits(x, "numeric")) {
      x <- ee$Number(x)
    }

    x$divide(e2)
  }
}

#' @keywords internal
`abs.ee.ee_number.Number` <- function(x) {
  if (!is.null(ee)) {
    if (inherits(x, "numeric")) {
      x <- ee$Number(x)
    }

    x$abs()
  }
}
#' @keywords internal
`acos.ee.ee_number.Number` <- function(x) {
  if (!is.null(ee)) {
    if (inherits(x, "numeric")) {
      x <- ee$Number(x)
    }

    x$acos()
  }
}
#' @keywords internal
`asin.ee.ee_number.Number` <- function(x) {
  if (!is.null(ee)) {
    if (inherits(x, "numeric")) {
      x <- ee$Number(x)
    }

    x$asin()
  }
}
#' @keywords internal
`atan.ee.ee_number.Number` <- function(x) {
  if (!is.null(ee)) {
    if (inherits(x, "numeric")) {
      x <- ee$Number(x)
    }

    x$atan()
  }
}

# Calculates the angle formed by the 2D vector x, y.
#
# @param x The number to calculate the `atan2` on.
#
# @return The resulting ee.Number with the results from `atan2`.
#
#' @keywords internal
`atan2.ee.ee_number.Number` <- function(x) {
  if (!is.null(ee)) {
    if (inherits(x, "numeric")) {
      x <- ee$Number(x)
    }

    x$atan2()
  }
}
#' @keywords internal
`floor.ee.ee_number.Number` <- function(x) {
  if (!is.null(ee)) {
    if (inherits(x, "numeric")) {
      x <- ee$Number(x)
    }

    x$floor()
  }
}
#' @keywords internal
`ceiling.ee.ee_number.Number` <- function(x) {
  if (!is.null(ee)) {
    if (inherits(x, "numeric")) {
      x <- ee$Number(x)
    }

    x$ceil()
  }
}
#' @keywords internal
`&.ee.ee_number.Number` <- function(x, e2) {
  if (!is.null(ee)) {
    if (inherits(x, "numeric")) {
      x <- ee$Number(x)
    }

    x$bitwiseAnd(e2)
  }
}
#' @keywords internal
`!.ee.ee_number.Number` <- function(x) {
  x$bitwiseNot()
}
#' @keywords internal
`|.ee.ee_number.Number` <- function(x, e2) {
  if (!is.null(ee)) {
    if (inherits(x, "numeric")) {
      x <- ee$Number(x)
    }

    x$bitwiseOr(e2)
  }
}
#' @keywords internal
`^.ee.ee_number.Number` <- function(x, e2) {
  if (!is.null(ee)) {
    if (inherits(x, "numeric")) {
      x <- ee$Number(x)
    }

    x$bitwiseXor(e2)
  }
}
#' @keywords internal
`cos.ee.ee_number.Number` <- function(x) {
  if (!is.null(ee)) {
    if (inherits(x, "numeric")) {
      x <- ee$Number(x)
    }

    x$cos()
  }
}
#' @keywords internal
`cosh.ee.ee_number.Number` <- function(x) {
  if (!is.null(ee)) {
    if (inherits(x, "numeric")) {
      x <- ee$Number(x)
    }

    x$cosh()
  }
}
#' @keywords internal
`digamma.ee.ee_number.Number` <- function(x) {
  if (!is.null(ee)) {
    if (inherits(x, "numeric")) {
      x <- ee$Number(x)
    }

    x$digamma()
  }
}
#' @keywords internal
`exp.ee.ee_number.Number` <- function(x) {
  if (!is.null(ee)) {
    if (inherits(x, "numeric")) {
      x <- ee$Number(x)
    }

    x$exp()
  }
}
#' @keywords internal
`as.numeric.ee.ee_number.Number` <- function(x, ...) {
  if (!is.null(ee)) {
    if (inherits(x, "numeric")) {
      x <- ee$Number(x)
    }

    x$double()
  }
}
#' @keywords internal
`gamma.ee.ee_number.Number` <- function(x) {
  if (!is.null(ee)) {
    if (inherits(x, "numeric")) {
      x <- ee$Number(x)
    }

    x$gamma()
  }
}
#' @keywords internal
`>.ee.ee_number.Number` <- function(x, e2) {
  if (!is.null(ee)) {
    if (inherits(x, "numeric")) {
      x <- ee$Number(x)
    }

    x$gt(e2)
  }
}
#' @keywords internal
`>=.ee.ee_number.Number` <- function(x, e2) {
  if (!is.null(ee)) {
    if (inherits(x, "numeric")) {
      x <- ee$Number(x)
    }

    x$gte(e2)
  }
}
#' @keywords internal
`<.ee.ee_number.Number` <- function(x, e2) {
  if (!is.null(ee)) {
    if (inherits(x, "numeric")) {
      x <- ee$Number(x)
    }

    x$lt(e2)
  }
}
#' @keywords internal
`<=.ee.ee_number.Number` <- function(x, e2) {
  if (!is.null(ee)) {
    if (inherits(x, "numeric")) {
      x <- ee$Number(x)
    }

    x$lte(e2)
  }
}
#' @keywords internal
`==.ee.ee_number.Number` <- function(x, e2) {
  if (!is.null(ee)) {
    if (inherits(x, "numeric")) {
      x <- ee$Number(x)
    }

    x$eq(e2)
  }
}
#' @keywords internal
`!=.ee.ee_number.Number` <- function(x, e2) {
  if (!is.null(ee)) {
    if (inherits(x, "numeric")) {
      x <- ee$Number(x)
    }

    x$neq(e2)
  }
}
#' @keywords internal
`as.integer.ee.ee_number.Number` <- function(x, ...) {
  if (!is.null(ee)) {
    if (inherits(x, "numeric")) {
      x <- ee$Number(x)
    }

    x$int()
  }
}
#' @keywords internal
`log.ee.ee_number.Number` <- function(x, base = NULL) {
  if (!is.null(ee)) {
    if (inherits(x, "numeric")) {
      x <- ee$Number(x)
    }

    x$log()
  }
}
#' @keywords internal
`log10.ee.ee_number.Number` <- function(x) {
  if (!is.null(ee)) {
    if (inherits(x, "numeric")) {
      x <- ee$Number(x)
    }

    x$log10()
  }
}
#' @keywords internal
`max.ee.ee_number.Number` <- function(x, ..., na.rm = FALSE) {
  args <- list(...)
  stopifnot("Expected exactly two arguments for max within GEE" = length(args) == 1)
  x <- args[[1]]
  if (!is.null(ee)) {
    if (inherits(x, "numeric")) {
      x <- ee$Number(x)
    }

    x$max(args[[2]])
  }
}
#' @keywords internal
`min.ee.ee_number.Number` <- function(x, ..., na.rm = FALSE) {
  args <- list(...)
  stopifnot("Expected exactly two arguments for min within GEE" = length(args) == 1)

  x$min(args[[1]])
}
#' @keywords internal
`%%.ee.ee_number.Number` <- function(x, e2) {
  if (!is.null(ee)) {
    if (inherits(x, "numeric")) {
      x <- ee$Number(x)
    }

    x$mod(e2)
  }
}
#' @keywords internal
`^.ee.ee_number.Number` <- function(x, e2) {
  if (!is.null(ee)) {
    if (inherits(x, "numeric")) {
      x <- ee$Number(x)
    }

    x$pow(e2)
  }
}

#' @keywords internal
`round.ee.ee_number.Number` <- function(x, digits = NULL) {
  if (!is.null(ee)) {
    if (inherits(x, "numeric")) {
      x <- ee$Number(x)
    }

    x$round()
  }
}
#' @keywords internal
`sin.ee.ee_number.Number` <- function(x) {
  if (!is.null(ee)) {
    if (inherits(x, "numeric")) {
      x <- ee$Number(x)
    }

    x$sin()
  }
}
#' @keywords internal
`sinh.ee.ee_number.Number` <- function(x) {
  if (!is.null(ee)) {
    if (inherits(x, "numeric")) {
      x <- ee$Number(x)
    }

    x$sinh()
  }
}
#' @keywords internal
`sqrt.ee.ee_number.Number` <- function(x) {
  if (!is.null(ee)) {
    if (inherits(x, "numeric")) {
      x <- ee$Number(x)
    }

    x$sqrt()
  }
}
#' @keywords internal
`tan.ee.ee_number.Number` <- function(x) {
  if (!is.null(ee)) {
    if (inherits(x, "numeric")) {
      x <- ee$Number(x)
    }

    x$tan()
  }
}
#' @keywords internal
`tanh.ee.ee_number.Number` <- function(x) {
  if (!is.null(ee)) {
    if (inherits(x, "numeric")) {
      x <- ee$Number(x)
    }

    x$tanh()
  }
}
#' @keywords internal
`trigamma.ee.ee_number.Number` <- function(x) {
  if (!is.null(ee)) {
    if (inherits(x, "numeric")) {
      x <- ee$Number(x)
    }

    x$trigamma()
  }
}
