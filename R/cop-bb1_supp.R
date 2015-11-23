## ----- Reflected -----

## Distribution

#' BB1 copula: supplemental functions
#'
#' Appended "u" means a horizontal flip (i.e. switching u to 1-u).
#' Appended "v" means a vertical flip (i.e. switching v to 1-v).
#' Appended "r" means reflection/survival copula
#'
#' @rdname bb1_supp
#' @export
pbb1r <- function(u, v, cpar) u + v - 1 + pbb1(1-u, 1-v, cpar)

## Density

#' @rdname bb1_supp
#' @export
dbb1r <- function(u, v, cpar) dbb1(1-u, 1-v, cpar)
#' @rdname bb1_supp
#' @export
logdbb1r <- function(u, v, cpar) logdbb1(1-u, 1-v, cpar)

## Random number generator

#' @rdname bb1_supp
#' @export
rbb1r <- function(n, cpar) 1 - rbb1(n, cpar)

## Conditional distribution

#' @rdname bb1_supp
#' @export
pcondbb1r <- function(v, u, cpar) 1 - pcondbb1(1-v, 1-u, cpar)
#' @rdname bb1_supp
#' @export
qcondbb1r <- function(p, u, cpar) 1 - qcondbb1(1-p, 1-u, cpar)

## ----- U-Flipped / V-axis reflection / Horizontal reflection -----

## Distribution

#' @rdname bb1_supp
#' @export
pbb1u <- function(u, v, cpar) v - pbb1(1-u, v, cpar)

## Density

#' @rdname bb1_supp
#' @export
dbb1u <- function(u, v, cpar) dbb1(1-u, v, cpar)
#' @rdname bb1_supp
#' @export
logdbb1u <- function(u, v, cpar) logdbb1(1-u, v, cpar)

## Random number generator

#' @rdname bb1_supp
#' @export
rbb1u <- function(n, cpar) {
    res <- rbb1(n, cpar)
    res[, 1] <- 1 - res[, 1]
    res
}

## Conditional distribution

#' @rdname bb1_supp
#' @export
pcondbb1u <- function(v, u, cpar) pcondbb1(v, 1-u, cpar)
#' @rdname bb1_supp
#' @export
qcondbb1u <- function(p, u, cpar) qcondbb1(p, 1-u, cpar)

## ----- V-Flipped / U-axis reflection / Vertical reflection -----

## Distribution

#' @rdname bb1_supp
#' @export
pbb1v <- function(u, v, cpar) u - pbb1(u, 1-v, cpar)

## Density

#' @rdname bb1_supp
#' @export
dbb1v <- function(u, v, cpar) dbb1(u, 1-v, cpar)
#' @rdname bb1_supp
#' @export
logdbb1v <- function(u, v, cpar) logdbb1(u, 1-v, cpar)

## Random number generator

#' @rdname bb1_supp
#' @export
rbb1v <- function(n, cpar) {
    res <- rbb1(n, cpar)
    res[, 2] <- 1 - res[, 2]
    res
}

## Conditional distribution

#' @rdname bb1_supp
#' @export
pcondbb1v <- function(v, u, cpar) 1 - pcondbb1(1-v, u, cpar)
#' @rdname bb1_supp
#' @export
qcondbb1v <- function(p, u, cpar) 1 - qcondbb1(1-p, u, cpar)

