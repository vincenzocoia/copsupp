## ----- U-Flipped / V-axis reflection / Horizontal reflection -----

## Distribution

#' MTCJ copula: supplemental functions
#'
#' Appended "u" means a horizontal flip (i.e. switching u to 1-u).
#' Appended "v" means a vertical flip (i.e. switching v to 1-v).
#' Appended "r" means reflection/survival copula
#'
#' @rdname mtcj_supp
#' @export
pmtcju <- function(u, v, cpar) v - CopulaModel::pmtcj(1-u, v, cpar)

## Density

#' @rdname mtcj_supp
#' @export
dmtcju <- function(u, v, cpar) CopulaModel::dmtcj(1-u, v, cpar)
#' @rdname mtcj_supp
#' @export
logdmtcju <- function(u, v, cpar) CopulaModel::logdmtcj(1-u, v, cpar)

## Random number generator

#' @rdname mtcj_supp
#' @export
rmtcju <- function(n, cpar) {
    res <- CopulaModel::rmtcj(n, cpar)
    res[, 1] <- 1 - res[, 1]
    res
}

## Conditional distribution

#' @rdname mtcj_supp
#' @export
pcondmtcju <- function(v, u, cpar) CopulaModel::pcondmtcj(v, 1-u, cpar)
#' @rdname mtcj_supp
#' @export
qcondmtcju <- function(p, u, cpar) CopulaModel::qcondmtcj(p, 1-u, cpar)

## ----- V-Flipped / U-axis reflection / Vertical reflection -----

## Distribution

#' @rdname mtcj_supp
#' @export
pmtcjv <- function(u, v, cpar) u - CopulaModel::pmtcj(u, 1-v, cpar)

## Density

#' @rdname mtcj_supp
#' @export
dmtcjv <- function(u, v, cpar) CopulaModel::dmtcj(u, 1-v, cpar)
#' @rdname mtcj_supp
#' @export
logdmtcjv <- function(u, v, cpar) CopulaModel::logdmtcj(u, 1-v, cpar)

## Random number generator

#' @rdname mtcj_supp
#' @export
rmtcjv <- function(n, cpar) {
    res <- CopulaModel::rmtcj(n, cpar)
    res[, 2] <- 1 - res[, 2]
    res
}

## Conditional distribution

#' @rdname mtcj_supp
#' @export
pcondmtcjv <- function(v, u, cpar) 1 - CopulaModel::pcondmtcj(1-v, u, cpar)
#' @rdname mtcj_supp
#' @export
qcondmtcjv <- function(p, u, cpar) 1 - CopulaModel::qcondmtcj(1-p, u, cpar)

