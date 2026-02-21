#' Creates an `Earth Engine` server number
#'
#' @param x the number to convert to Earth Engine's number.
#'
#' @return The Earth Engine number
#'
#' @seealso https://developers.google.com/earth-engine/apidocs/ee-number
#' @export
ee_number <- function(x) {
  UseMethod("ee_number")
}

#' @export
"ee_number.numeric" <- function(x) if (!is.null(ee)) ee$Number(x)

#' @export
`+.ee.ee_number.Number` <- function(x, e2) {
  if (!is.null(ee)) {
    if (inherits(x, "numeric")) {
      x <- ee$Number(x)
    }

    x$add(e2)
  }
}

#' @export
`-.ee.ee_number.Number` <- function(x, e2) {
  if (!is.null(ee)) {
    if (inherits(x, "numeric")) {
      x <- ee$Number(x)
    }

    x$subtract(e2)
  }
}

#' @export
`*.ee.ee_number.Number` <- function(x, e2) {
  if (!is.null(ee)) {
    if (inherits(x, "numeric")) {
      x <- ee$Number(x)
    }

    x$multiply(e2)
  }
}

#' @export
`/.ee.ee_number.Number` <- function(x, e2) {
  if (!is.null(ee)) {
    if (inherits(x, "numeric")) {
      x <- ee$Number(x)
    }

    x$divide(e2)
  }
}

#' @export
`abs.ee.ee_number.Number` <- function(x) {
  if (!is.null(ee)) {
    if (inherits(x, "numeric")) {
      x <- ee$Number(x)
    }

    x$abs()
  }
}
#' @export
`acos.ee.ee_number.Number` <- function(x) {
  if (!is.null(ee)) {
    if (inherits(x, "numeric")) {
      x <- ee$Number(x)
    }

    x$acos()
  }
}
#' @export
`asin.ee.ee_number.Number` <- function(x) {
  if (!is.null(ee)) {
    if (inherits(x, "numeric")) {
      x <- ee$Number(x)
    }

    x$asin()
  }
}
#' @export
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
#' @export
`atan2.ee.ee_number.Number` <- function(x) {
  if (!is.null(ee)) {
    if (inherits(x, "numeric")) {
      x <- ee$Number(x)
    }

    x$atan2()
  }
}
#' @export
`floor.ee.ee_number.Number` <- function(x) {
  if (!is.null(ee)) {
    if (inherits(x, "numeric")) {
      x <- ee$Number(x)
    }

    x$floor()
  }
}
#' @export
`ceiling.ee.ee_number.Number` <- function(x) {
  if (!is.null(ee)) {
    if (inherits(x, "numeric")) {
      x <- ee$Number(x)
    }

    x$ceil()
  }
}
#' @export
`&.ee.ee_number.Number` <- function(x, e2) {
  if (!is.null(ee)) {
    if (inherits(x, "numeric")) {
      x <- ee$Number(x)
    }

    x$bitwiseAnd(e2)
  }
}
#' @export
`!.ee.ee_number.Number` <- function(x) {
  x$bitwiseNot()
}
#' @export
`|.ee.ee_number.Number` <- function(x, e2) {
  if (!is.null(ee)) {
    if (inherits(x, "numeric")) {
      x <- ee$Number(x)
    }

    x$bitwiseOr(e2)
  }
}
#' @export
`^.ee.ee_number.Number` <- function(x, e2) {
  if (!is.null(ee)) {
    if (inherits(x, "numeric")) {
      x <- ee$Number(x)
    }

    x$bitwiseXor(e2)
  }
}
#' @export
`cos.ee.ee_number.Number` <- function(x) {
  if (!is.null(ee)) {
    if (inherits(x, "numeric")) {
      x <- ee$Number(x)
    }

    x$cos()
  }
}
#' @export
`cosh.ee.ee_number.Number` <- function(x) {
  if (!is.null(ee)) {
    if (inherits(x, "numeric")) {
      x <- ee$Number(x)
    }

    x$cosh()
  }
}
#' @export
`digamma.ee.ee_number.Number` <- function(x) {
  if (!is.null(ee)) {
    if (inherits(x, "numeric")) {
      x <- ee$Number(x)
    }

    x$digamma()
  }
}
#' @export
`exp.ee.ee_number.Number` <- function(x) {
  if (!is.null(ee)) {
    if (inherits(x, "numeric")) {
      x <- ee$Number(x)
    }

    x$exp()
  }
}
#' @export
`as.numeric.ee.ee_number.Number` <- function(x, ...) {
  if (!is.null(ee)) {
    if (inherits(x, "numeric")) {
      x <- ee$Number(x)
    }

    x$double()
  }
}
#' @export
`gamma.ee.ee_number.Number` <- function(x) {
  if (!is.null(ee)) {
    if (inherits(x, "numeric")) {
      x <- ee$Number(x)
    }

    x$gamma()
  }
}
#' @export
`>.ee.ee_number.Number` <- function(x, e2) {
  if (!is.null(ee)) {
    if (inherits(x, "numeric")) {
      x <- ee$Number(x)
    }

    x$gt(e2)
  }
}
#' @export
`>=.ee.ee_number.Number` <- function(x, e2) {
  if (!is.null(ee)) {
    if (inherits(x, "numeric")) {
      x <- ee$Number(x)
    }

    x$gte(e2)
  }
}
#' @export
`<.ee.ee_number.Number` <- function(x, e2) {
  if (!is.null(ee)) {
    if (inherits(x, "numeric")) {
      x <- ee$Number(x)
    }

    x$lt(e2)
  }
}
#' @export
`<=.ee.ee_number.Number` <- function(x, e2) {
  if (!is.null(ee)) {
    if (inherits(x, "numeric")) {
      x <- ee$Number(x)
    }

    x$lte(e2)
  }
}
#' @export
`==.ee.ee_number.Number` <- function(x, e2) {
  if (!is.null(ee)) {
    if (inherits(x, "numeric")) {
      x <- ee$Number(x)
    }

    x$eq(e2)
  }
}
#' @export
`!=.ee.ee_number.Number` <- function(x, e2) {
  if (!is.null(ee)) {
    if (inherits(x, "numeric")) {
      x <- ee$Number(x)
    }

    x$neq(e2)
  }
}
#' @export
`as.integer.ee.ee_number.Number` <- function(x, ...) {
  if (!is.null(ee)) {
    if (inherits(x, "numeric")) {
      x <- ee$Number(x)
    }

    x$int()
  }
}
#' @export
`log.ee.ee_number.Number` <- function(x, base = NULL) {
  if (!is.null(ee)) {
    if (inherits(x, "numeric")) {
      x <- ee$Number(x)
    }

    x$log()
  }
}
#' @export
`log10.ee.ee_number.Number` <- function(x) {
  if (!is.null(ee)) {
    if (inherits(x, "numeric")) {
      x <- ee$Number(x)
    }

    x$log10()
  }
}
#' @export
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
#' @export
`min.ee.ee_number.Number` <- function(x, ..., na.rm = FALSE) {
  args <- list(...)
  stopifnot("Expected exactly two arguments for min within GEE" = length(args) == 1)

  x$min(args[[1]])
}
#' @export
`%%.ee.ee_number.Number` <- function(x, e2) {
  if (!is.null(ee)) {
    if (inherits(x, "numeric")) {
      x <- ee$Number(x)
    }

    x$mod(e2)
  }
}
#' @export
`^.ee.ee_number.Number` <- function(x, e2) {
  if (!is.null(ee)) {
    if (inherits(x, "numeric")) {
      x <- ee$Number(x)
    }

    x$pow(e2)
  }
}

#' @export
`round.ee.ee_number.Number` <- function(x, digits = NULL) {
  if (!is.null(ee)) {
    if (inherits(x, "numeric")) {
      x <- ee$Number(x)
    }

    x$round()
  }
}
#' @export
`sin.ee.ee_number.Number` <- function(x) {
  if (!is.null(ee)) {
    if (inherits(x, "numeric")) {
      x <- ee$Number(x)
    }

    x$sin()
  }
}
#' @export
`sinh.ee.ee_number.Number` <- function(x) {
  if (!is.null(ee)) {
    if (inherits(x, "numeric")) {
      x <- ee$Number(x)
    }

    x$sinh()
  }
}
#' @export
`sqrt.ee.ee_number.Number` <- function(x) {
  if (!is.null(ee)) {
    if (inherits(x, "numeric")) {
      x <- ee$Number(x)
    }

    x$sqrt()
  }
}
#' @export
`tan.ee.ee_number.Number` <- function(x) {
  if (!is.null(ee)) {
    if (inherits(x, "numeric")) {
      x <- ee$Number(x)
    }

    x$tan()
  }
}
#' @export
`tanh.ee.ee_number.Number` <- function(x) {
  if (!is.null(ee)) {
    if (inherits(x, "numeric")) {
      x <- ee$Number(x)
    }

    x$tanh()
  }
}
#' @export
`trigamma.ee.ee_number.Number` <- function(x) {
  if (!is.null(ee)) {
    if (inherits(x, "numeric")) {
      x <- ee$Number(x)
    }

    x$trigamma()
  }
}
