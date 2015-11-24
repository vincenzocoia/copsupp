## ----- Reflected -----

## Distribution

#' BB7 copula: supplemental functions
#'
#' Appended "u" means a horizontal flip (i.e. switching u to 1-u).
#' Appended "v" means a vertical flip (i.e. switching v to 1-v).
#' Appended "r" means reflection/survival copula
#'
#' @rdname bb7_supp
#' @export
pbb7r <- function(u, v, cpar) u + v - 1 + CopulaModel::pbb7(1-u, 1-v, cpar)

## Density

#' @rdname bb7_supp
#' @export
dbb7r <- function(u, v, cpar) CopulaModel::dbb7(1-u, 1-v, cpar)
#' @rdname bb7_supp
#' @export
logdbb7r <- function(u, v, cpar) CopulaModel::logdbb7(1-u, 1-v, cpar)

## Random number generator

#' @rdname bb7_supp
#' @export
rbb7r <- function(n, cpar) 1 - CopulaModel::rbb7(n, cpar)

## Conditional distribution

#' @rdname bb7_supp
#' @export
pcondbb7r <- function(v, u, cpar) 1 - CopulaModel::pcondbb7(1-v, 1-u, cpar)
#' @rdname bb7_supp
#' @export
qcondbb7r <- function(p, u, cpar) 1 - CopulaModel::qcondbb7(1-p, 1-u, cpar)

## ----- U-Flipped / V-axis reflection / Horizontal reflection -----

## Distribution

#' @rdname bb7_supp
#' @export
pbb7u <- function(u, v, cpar) v - CopulaModel::pbb7(1-u, v, cpar)

## Density

#' @rdname bb7_supp
#' @export
dbb7u <- function(u, v, cpar) CopulaModel::dbb7(1-u, v, cpar)
#' @rdname bb7_supp
#' @export
logdbb7u <- function(u, v, cpar) CopulaModel::logdbb7(1-u, v, cpar)

## Random number generator

#' @rdname bb7_supp
#' @export
rbb7u <- function(n, cpar) {
    res <- CopulaModel::rbb7(n, cpar)
    res[, 1] <- 1 - res[, 1]
    res
}

## Conditional distribution

#' @rdname bb7_supp
#' @export
pcondbb7u <- function(v, u, cpar) CopulaModel::pcondbb7(v, 1-u, cpar)
#' @rdname bb7_supp
#' @export
qcondbb7u <- function(p, u, cpar) CopulaModel::qcondbb7(p, 1-u, cpar)

## ----- V-Flipped / U-axis reflection / Vertical reflection -----

## Distribution

#' @rdname bb7_supp
#' @export
pbb7v <- function(u, v, cpar) u - CopulaModel::pbb7(u, 1-v, cpar)

## Density

#' @rdname bb7_supp
#' @export
dbb7v <- function(u, v, cpar) CopulaModel::dbb7(u, 1-v, cpar)
#' @rdname bb7_supp
#' @export
logdbb7v <- function(u, v, cpar) CopulaModel::logdbb7(u, 1-v, cpar)

## Random number generator

#' @rdname bb7_supp
#' @export
rbb7v <- function(n, cpar) {
    res <- CopulaModel::rbb7(n, cpar)
    res[, 2] <- 1 - res[, 2]
    res
}

## Conditional distribution

#' @rdname bb7_supp
#' @export
pcondbb7v <- function(v, u, cpar) 1 - CopulaModel::pcondbb7(1-v, u, cpar)
#' @rdname bb7_supp
#' @export
qcondbb7v <- function(p, u, cpar) 1 - CopulaModel::qcondbb7(1-p, u, cpar)

