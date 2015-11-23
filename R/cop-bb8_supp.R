## ----- Reflected -----

## Distribution

#' BB8 copula: supplemental functions
#'
#' Appended "u" means a horizontal flip (i.e. switching u to 1-u).
#' Appended "v" means a vertical flip (i.e. switching v to 1-v).
#' Appended "r" means reflection/survival copula
#'
#' @rdname bb8_supp
#' @export
pbb8r <- function(u, v, cpar) u + v - 1 + pbb8(1-u, 1-v, cpar)

## Density

#' @rdname bb8_supp
#' @export
dbb8r <- function(u, v, cpar) dbb8(1-u, 1-v, cpar)
#' @rdname bb8_supp
#' @export
logdbb8r <- function(u, v, cpar) logdbb8(1-u, 1-v, cpar)

## Random number generator

#' @rdname bb8_supp
#' @export
rbb8r <- function(n, cpar) 1 - rbb8(n, cpar)

## Conditional distribution

#' @rdname bb8_supp
#' @export
pcondbb8r <- function(v, u, cpar) 1 - pcondbb8(1-v, 1-u, cpar)
#' @rdname bb8_supp
#' @export
qcondbb8r <- function(p, u, cpar) 1 - qcondbb8(1-p, 1-u, cpar)

## ----- U-Flipped / V-axis reflection / Horizontal reflection -----

## Distribution

#' @rdname bb8_supp
#' @export
pbb8u <- function(u, v, cpar) v - pbb8(1-u, v, cpar)

## Density

#' @rdname bb8_supp
#' @export
dbb8u <- function(u, v, cpar) dbb8(1-u, v, cpar)
#' @rdname bb8_supp
#' @export
logdbb8u <- function(u, v, cpar) logdbb8(1-u, v, cpar)

## Random number generator

#' @rdname bb8_supp
#' @export
rbb8u <- function(n, cpar) {
    res <- rbb8(n, cpar)
    res[, 1] <- 1 - res[, 1]
    res
}

## Conditional distribution

#' @rdname bb8_supp
#' @export
pcondbb8u <- function(v, u, cpar) pcondbb8(v, 1-u, cpar)
#' @rdname bb8_supp
#' @export
qcondbb8u <- function(p, u, cpar) qcondbb8(p, 1-u, cpar)

## ----- V-Flipped / U-axis reflection / Vertical reflection -----

## Distribution

#' @rdname bb8_supp
#' @export
pbb8v <- function(u, v, cpar) u - pbb8(u, 1-v, cpar)

## Density

#' @rdname bb8_supp
#' @export
dbb8v <- function(u, v, cpar) dbb8(u, 1-v, cpar)
#' @rdname bb8_supp
#' @export
logdbb8v <- function(u, v, cpar) logdbb8(u, 1-v, cpar)

## Random number generator

#' @rdname bb8_supp
#' @export
rbb8v <- function(n, cpar) {
    res <- rbb8(n, cpar)
    res[, 2] <- 1 - res[, 2]
    res
}

## Conditional distribution

#' @rdname bb8_supp
#' @export
pcondbb8v <- function(v, u, cpar) 1 - pcondbb8(1-v, u, cpar)
#' @rdname bb8_supp
#' @export
qcondbb8v <- function(p, u, cpar) 1 - qcondbb8(1-p, u, cpar)

