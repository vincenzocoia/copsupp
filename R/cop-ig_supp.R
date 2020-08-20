

## ----- U-Flipped / V-axis reflection / Horizontal reflection -----

## Distribution

#' IG copula: supplemental functions
#'
#' Appended "u" means a horizontal flip (i.e. switching u to 1-u).
#' Appended "v" means a vertical flip (i.e. switching v to 1-v).
#' Appended "r" means reflection/survival copula
#'
#' @rdname igcop_supp
#' @export
pigcopu <- function(u, v, cpar) v - pigcop(1-u, v, cpar)

## Density

#' @rdname igcop_supp
#' @export
digcopu <- function(u, v, cpar) digcop(1-u, v, cpar)
#' @rdname igcop_supp
#' @export
logdigcopu <- function(u, v, cpar) logdigcop(1-u, v, cpar)

## Random number generator

#' @rdname igcop_supp
#' @export
rigcopu <- function(n, cpar) {
    res <- rigcop(n, cpar)
    res[, 1] <- 1 - res[, 1]
    res
}

## Conditional distribution

#' @rdname igcop_supp
#' @export
pcondigcopu <- function(v, u, cpar) pcondigcop(v, 1-u, cpar)
#' @rdname igcop_supp
#' @export
qcondigcopu <- function(p, u, cpar) qcondigcop(p, 1-u, cpar)

## ----- V-Flipped / U-axis reflection / Vertical reflection -----

## Distribution

#' @rdname igcop_supp
#' @export
pigcopv <- function(u, v, cpar) u - pigcop(u, 1-v, cpar)

## Density

#' @rdname igcop_supp
#' @export
digcopv <- function(u, v, cpar) digcop(u, 1-v, cpar)
#' @rdname igcop_supp
#' @export
logdigcopv <- function(u, v, cpar) logdigcop(u, 1-v, cpar)

## Conditional distribution

#' @rdname igcop_supp
#' @export
pcondigcopv <- function(v, u, cpar) 1 - pcondigcop(1-v, u, cpar)
#' @rdname igcop_supp
#' @export
qcondigcopv <- function(p, u, cpar) 1 - qcondigcop(1-p, u, cpar)

## ----- Reflection ------

#' @rdname igcop_supp
#' @export
pcondigcopr <- function(v, u, cpar) 1 - pcondigcop(1 - v, 1 - u, cpar)

#' @rdname igcop_supp
#' @export
qcondigcopr <- function(p, u, cpar) 1 - qcondigcop(1 - p, 1 - u, cpar)

#' @rdname igcop_supp
#' @export
pcondigcopr12 <- function(u, v, cpar) 1 - pcondigcop12(1 - u, 1 - v, cpar)

#' @rdname igcop_supp
#' @export
qcondigcopr12 <- function(p, v, cpar) 1 - qcondigcop12(1 - p, 1 - v, cpar)

#' @rdname igcop_supp
#' @export
digcopr <- function(u, v, cpar) digcop(1 - u, 1 - v, cpar)

#' @rdname igcop_supp
#' @export
pigcopr <- function(u, v, cpar) u + v - 1 + pigcop(1 - u, 1 - v, cpar)

#' @rdname igcop_supp
#' @export
rigcopr <- function(n, cpar) 1 - rigcop(n, cpar)
