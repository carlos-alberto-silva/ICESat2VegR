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
`+.ee.ee_number.Number` <- function(e1, e2) {
  if (!is.null(ee)) {
    if (inherits(e1, "numeric")) {
      e1 <- ee$Number(e1)
    }

    e1$add(e2)
  }
}

#' @export
`-.ee.ee_number.Number` <- function(e1, e2) {
  if (!is.null(ee)) {
    if (inherits(e1, "numeric")) {
      e1 <- ee$Number(e1)
    }

    e1$subtract(e2)
  }
}

#' @export
`*.ee.ee_number.Number` <- function(e1, e2) {
  if (!is.null(ee)) {
    if (inherits(e1, "numeric")) {
      e1 <- ee$Number(e1)
    }

    e1$multiply(e2)
  }
}

#' @export
`/.ee.ee_number.Number` <- function(e1, e2) {
  if (!is.null(ee)) {
    if (inherits(e1, "numeric")) {
      e1 <- ee$Number(e1)
    }

    e1$divide(e2)
  }
}

#' @export
`abs.ee.ee_number.Number` <- function(e1) {
  if (!is.null(ee)) {
    if (inherits(e1, "numeric")) {
      e1 <- ee$Number(e1)
    }

    e1$abs()
  }
}
#' @export
`acos.ee.ee_number.Number` <- function(e1) {
  if (!is.null(ee)) {
    if (inherits(e1, "numeric")) {
      e1 <- ee$Number(e1)
    }

    e1$acos()
  }
}
#' @export
`asin.ee.ee_number.Number` <- function(e1) {
  if (!is.null(ee)) {
    if (inherits(e1, "numeric")) {
      e1 <- ee$Number(e1)
    }

    e1$asin()
  }
}
#' @export
`atan.ee.ee_number.Number` <- function(e1) {
  if (!is.null(ee)) {
    if (inherits(e1, "numeric")) {
      e1 <- ee$Number(e1)
    }

    e1$atan()
  }
}

#' Calculates the angle formed by the 2D vector x, y.
#'
#' @param e1 The number to calculate the `atan2` on.
#'
#' @return The resulting ee.Number with the results from `atan2`.
#'
#' @export
`atan2.ee.ee_number.Number` <- function(e1) {
  if (!is.null(ee)) {
    if (inherits(e1, "numeric")) {
      e1 <- ee$Number(e1)
    }

    e1$atan2()
  }
}
#' @export
`floor.ee.ee_number.Number` <- function(e1) {
  if (!is.null(ee)) {
    if (inherits(e1, "numeric")) {
      e1 <- ee$Number(e1)
    }

    e1$floor()
  }
}
#' @export
`ceiling.ee.ee_number.Number` <- function(e1) {
  if (!is.null(ee)) {
    if (inherits(e1, "numeric")) {
      e1 <- ee$Number(e1)
    }

    e1$ceil()
  }
}
#' @export
`&.ee.ee_number.Number` <- function(e1, e2) {
  if (!is.null(ee)) {
    if (inherits(e1, "numeric")) {
      e1 <- ee$Number(e1)
    }

    e1$bitwiseAnd(e2)
  }
}
#' @export
`!.ee.ee_number.Number` <- function(e1, e2) {
  if (!is.null(ee)) {
    if (inherits(e1, "numeric")) {
      e1 <- ee$Number(e1)
    }

    e1$bitwiseNot(e2)
  }
}
#' @export
`|.ee.ee_number.Number` <- function(e1, e2) {
  if (!is.null(ee)) {
    if (inherits(e1, "numeric")) {
      e1 <- ee$Number(e1)
    }

    e1$bitwiseOr(e2)
  }
}
#' @export
`^.ee.ee_number.Number` <- function(e1, e2) {
  if (!is.null(ee)) {
    if (inherits(e1, "numeric")) {
      e1 <- ee$Number(e1)
    }

    e1$bitwiseXor(e2)
  }
}
#' @export
`cos.ee.ee_number.Number` <- function(e1) {
  if (!is.null(ee)) {
    if (inherits(e1, "numeric")) {
      e1 <- ee$Number(e1)
    }

    e1$cos()
  }
}
#' @export
`cosh.ee.ee_number.Number` <- function(e1) {
  if (!is.null(ee)) {
    if (inherits(e1, "numeric")) {
      e1 <- ee$Number(e1)
    }

    e1$cosh()
  }
}
#' @export
`digamma.ee.ee_number.Number` <- function(e1) {
  if (!is.null(ee)) {
    if (inherits(e1, "numeric")) {
      e1 <- ee$Number(e1)
    }

    e1$digamma()
  }
}
#' @export
`exp.ee.ee_number.Number` <- function(e1) {
  if (!is.null(ee)) {
    if (inherits(e1, "numeric")) {
      e1 <- ee$Number(e1)
    }

    e1$exp()
  }
}
#' @export
`as.numeric.ee.ee_number.Number` <- function(e1) {
  if (!is.null(ee)) {
    if (inherits(e1, "numeric")) {
      e1 <- ee$Number(e1)
    }

    e1$double()
  }
}
#' @export
`gamma.ee.ee_number.Number` <- function(e1) {
  if (!is.null(ee)) {
    if (inherits(e1, "numeric")) {
      e1 <- ee$Number(e1)
    }

    e1$gamma()
  }
}
#' @export
`>.ee.ee_number.Number` <- function(e1, e2) {
  if (!is.null(ee)) {
    if (inherits(e1, "numeric")) {
      e1 <- ee$Number(e1)
    }

    e1$gt(e2)
  }
}
#' @export
`>=.ee.ee_number.Number` <- function(e1, e2) {
  if (!is.null(ee)) {
    if (inherits(e1, "numeric")) {
      e1 <- ee$Number(e1)
    }

    e1$gte(e2)
  }
}
#' @export
`<.ee.ee_number.Number` <- function(e1, e2) {
  if (!is.null(ee)) {
    if (inherits(e1, "numeric")) {
      e1 <- ee$Number(e1)
    }

    e1$lt(e2)
  }
}
#' @export
`<=.ee.ee_number.Number` <- function(e1, e2) {
  if (!is.null(ee)) {
    if (inherits(e1, "numeric")) {
      e1 <- ee$Number(e1)
    }

    e1$lte(e2)
  }
}
#' @export
`==.ee.ee_number.Number` <- function(e1, e2) {
  if (!is.null(ee)) {
    if (inherits(e1, "numeric")) {
      e1 <- ee$Number(e1)
    }

    e1$eq(e2)
  }
}
#' @export
`!=.ee.ee_number.Number` <- function(e1, e2) {
  if (!is.null(ee)) {
    if (inherits(e1, "numeric")) {
      e1 <- ee$Number(e1)
    }

    e1$neq(e2)
  }
}
#' @export
`as.integer.ee.ee_number.Number` <- function(e1) {
  if (!is.null(ee)) {
    if (inherits(e1, "numeric")) {
      e1 <- ee$Number(e1)
    }

    e1$int()
  }
}
#' @export
`log.ee.ee_number.Number` <- function(e1) {
  if (!is.null(ee)) {
    if (inherits(e1, "numeric")) {
      e1 <- ee$Number(e1)
    }

    e1$log()
  }
}
#' @export
`log10.ee.ee_number.Number` <- function(e1) {
  if (!is.null(ee)) {
    if (inherits(e1, "numeric")) {
      e1 <- ee$Number(e1)
    }

    e1$log10()
  }
}
#' @export
`max.ee.ee_number.Number` <- function(e1, e2, na.rm = FALSE) {
  if (!is.null(ee)) {
    if (inherits(e1, "numeric")) {
      e1 <- ee$Number(e1)
    }

    e1$max(e2)
  }
}
#' @export
`min.ee.ee_number.Number` <- function(e1, e2, na.rm = FALSE) {
  if (!is.null(ee)) {
    if (inherits(e1, "numeric")) {
      e1 <- ee$Number(e1)
    }

    e1$min(e2)
  }
}
#' @export
`%%.ee.ee_number.Number` <- function(e1, e2) {
  if (!is.null(ee)) {
    if (inherits(e1, "numeric")) {
      e1 <- ee$Number(e1)
    }

    e1$mod(e2)
  }
}
#' @export
`^.ee.ee_number.Number` <- function(e1, e2) {
  if (!is.null(ee)) {
    if (inherits(e1, "numeric")) {
      e1 <- ee$Number(e1)
    }

    e1$pow(e2)
  }
}
#' @export
`round.ee.ee_number.Number` <- function(e1) {
  if (!is.null(ee)) {
    if (inherits(e1, "numeric")) {
      e1 <- ee$Number(e1)
    }

    e1$round()
  }
}
#' @export
`sin.ee.ee_number.Number` <- function(e1) {
  if (!is.null(ee)) {
    if (inherits(e1, "numeric")) {
      e1 <- ee$Number(e1)
    }

    e1$sin()
  }
}
#' @export
`sinh.ee.ee_number.Number` <- function(e1) {
  if (!is.null(ee)) {
    if (inherits(e1, "numeric")) {
      e1 <- ee$Number(e1)
    }

    e1$sinh()
  }
}
#' @export
`sqrt.ee.ee_number.Number` <- function(e1) {
  if (!is.null(ee)) {
    if (inherits(e1, "numeric")) {
      e1 <- ee$Number(e1)
    }

    e1$sqrt()
  }
}
#' @export
`tan.ee.ee_number.Number` <- function(e1) {
  if (!is.null(ee)) {
    if (inherits(e1, "numeric")) {
      e1 <- ee$Number(e1)
    }

    e1$tan()
  }
}
#' @export
`tanh.ee.ee_number.Number` <- function(e1) {
  if (!is.null(ee)) {
    if (inherits(e1, "numeric")) {
      e1 <- ee$Number(e1)
    }

    e1$tanh()
  }
}
#' @export
`trigamma.ee.ee_number.Number` <- function(e1) {
  if (!is.null(ee)) {
    if (inherits(e1, "numeric")) {
      e1 <- ee$Number(e1)
    }

    e1$trigamma()
  }
}
