## ----- Original ------

#' BB6 copula: supplemental functions
#'
#' Appended "u" means a horizontal flip (i.e. switching u to 1-u).
#' Appended "v" means a vertical flip (i.e. switching v to 1-v).
#' Appended "r" means reflection/survival copula
#'
#' @rdname bb6_supp
#' @export
logdbb6 <- function(u, v, cpar) log(CopulaModel::dbb6(u, v, cpar))

## ----- Reflected -----

## Distribution

#' @rdname bb6_supp
#' @export
pbb6r <- function(u, v, cpar) u + v - 1 + CopulaModel::pbb6(1-u, 1-v, cpar)

## Density

#' @rdname bb6_supp
#' @export
dbb6r <- function(u, v, cpar) CopulaModel::dbb6(1-u, 1-v, cpar)
#' @rdname bb6_supp
#' @export
logdbb6r <- function(u, v, cpar) CopulaModel::logdbb6(1-u, 1-v, cpar)

## Random number generator

#' @rdname bb6_supp
#' @export
rbb6r <- function(n, cpar) 1 - CopulaModel::rbb6(n, cpar)

## Conditional distribution

#' @rdname bb6_supp
#' @export
pcondbb6r <- function(v, u, cpar) 1 - CopulaModel::pcondbb6(1-v, 1-u, cpar)
#' @rdname bb6_supp
#' @export
qcondbb6r <- function(p, u, cpar) 1 - CopulaModel::qcondbb6(1-p, 1-u, cpar)

## ----- U-Flipped / V-axis reflection / Horizontal reflection -----

## Distribution

#' @rdname bb6_supp
#' @export
pbb6u <- function(u, v, cpar) v - CopulaModel::pbb6(1-u, v, cpar)

## Density

#' @rdname bb6_supp
#' @export
dbb6u <- function(u, v, cpar) CopulaModel::dbb6(1-u, v, cpar)
#' @rdname bb6_supp
#' @export
logdbb6u <- function(u, v, cpar) CopulaModel::logdbb6(1-u, v, cpar)

## Random number generator

#' @rdname bb6_supp
#' @export
rbb6u <- function(n, cpar) {
    res <- CopulaModel::rbb6(n, cpar)
    res[, 1] <- 1 - res[, 1]
    res
}

## Conditional distribution

#' @rdname bb6_supp
#' @export
pcondbb6u <- function(v, u, cpar) CopulaModel::pcondbb6(v, 1-u, cpar)
#' @rdname bb6_supp
#' @export
qcondbb6u <- function(p, u, cpar) CopulaModel::qcondbb6(p, 1-u, cpar)

## ----- V-Flipped / U-axis reflection / Vertical reflection -----

## Distribution

#' @rdname bb6_supp
#' @export
pbb6v <- function(u, v, cpar) u - CopulaModel::pbb6(u, 1-v, cpar)

## Density

#' @rdname bb6_supp
#' @export
dbb6v <- function(u, v, cpar) CopulaModel::dbb6(u, 1-v, cpar)
#' @rdname bb6_supp
#' @export
logdbb6v <- function(u, v, cpar) CopulaModel::logdbb6(u, 1-v, cpar)

## Random number generator

#' @rdname bb6_supp
#' @export
rbb6v <- function(n, cpar) {
    res <- CopulaModel::rbb6(n, cpar)
    res[, 2] <- 1 - res[, 2]
    res
}

## Conditional distribution

#' @rdname bb6_supp
#' @export
pcondbb6v <- function(v, u, cpar) 1 - CopulaModel::pcondbb6(1-v, u, cpar)
#' @rdname bb6_supp
#' @export
qcondbb6v <- function(p, u, cpar) 1 - CopulaModel::qcondbb6(1-p, u, cpar)

