## ----- Reflected -----

## Distribution

#' Joe copula: supplemental functions
#'
#' Appended "u" means a horizontal flip (i.e. switching u to 1-u).
#' Appended "v" means a vertical flip (i.e. switching v to 1-v).
#' Appended "r" means reflection/survival copula
#'
#' @rdname joe_supp
#' @export
pjoer <- function(u, v, cpar) u + v - 1 + CopulaModel::pjoe(1-u, 1-v, cpar)

## Density

#' @rdname joe_supp
#' @export
djoer <- function(u, v, cpar) CopulaModel::djoe(1-u, 1-v, cpar)
#' @rdname joe_supp
#' @export
logdjoer <- function(u, v, cpar) CopulaModel::logdjoe(1-u, 1-v, cpar)

## Random number generator

#' @rdname joe_supp
#' @export
rjoer <- function(n, cpar) 1 - CopulaModel::rjoe(n, cpar)

## Conditional distribution

#' @rdname joe_supp
#' @export
pcondjoer <- function(v, u, cpar) 1 - CopulaModel::pcondjoe(1-v, 1-u, cpar)
#' @rdname joe_supp
#' @export
qcondjoer <- function(p, u, cpar) 1 - CopulaModel::qcondjoe(1-p, 1-u, cpar)

## ----- U-Flipped / V-axis reflection / Horizontal reflection -----

## Distribution

#' @rdname joe_supp
#' @export
pjoeu <- function(u, v, cpar) v - CopulaModel::pjoe(1-u, v, cpar)

## Density

#' @rdname joe_supp
#' @export
djoeu <- function(u, v, cpar) CopulaModel::djoe(1-u, v, cpar)
#' @rdname joe_supp
#' @export
logdjoeu <- function(u, v, cpar) CopulaModel::logdjoe(1-u, v, cpar)

## Random number generator

#' @rdname joe_supp
#' @export
rjoeu <- function(n, cpar) {
    res <- CopulaModel::rjoe(n, cpar)
    res[, 1] <- 1 - res[, 1]
    res
}

## Conditional distribution

#' @rdname joe_supp
#' @export
pcondjoeu <- function(v, u, cpar) CopulaModel::pcondjoe(v, 1-u, cpar)
#' @rdname joe_supp
#' @export
qcondjoeu <- function(p, u, cpar) CopulaModel::qcondjoe(p, 1-u, cpar)

## ----- V-Flipped / U-axis reflection / Vertical reflection -----

## Distribution

#' @rdname joe_supp
#' @export
pjoev <- function(u, v, cpar) u - CopulaModel::pjoe(u, 1-v, cpar)

## Density

#' @rdname joe_supp
#' @export
djoev <- function(u, v, cpar) CopulaModel::djoe(u, 1-v, cpar)
#' @rdname joe_supp
#' @export
logdjoev <- function(u, v, cpar) CopulaModel::logdjoe(u, 1-v, cpar)

## Random number generator

#' @rdname joe_supp
#' @export
rjoev <- function(n, cpar) {
    res <- CopulaModel::rjoe(n, cpar)
    res[, 2] <- 1 - res[, 2]
    res
}

## Conditional distribution

#' @rdname joe_supp
#' @export
pcondjoev <- function(v, u, cpar) 1 - CopulaModel::pcondjoe(1-v, u, cpar)
#' @rdname joe_supp
#' @export
qcondjoev <- function(p, u, cpar) 1 - CopulaModel::qcondjoe(1-p, u, cpar)

