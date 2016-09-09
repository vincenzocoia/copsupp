## ----- U-Flipped / V-axis reflection / Horizontal reflection -----

## Distribution

#' Bivariate Skew Normal Copula: supplement
#'
#' @rdname bskewncop_supp
#' @export
pbskewncopu <- function(u, v, cpar) v - pbskewncop(1-u, v, cpar)

## Density

#' @rdname bskewncop_supp
#' @export
dbskewncopu <- function(u, v, cpar) dbskewncop(1-u, v, cpar)
#' @rdname bskewncop_supp
#' @export
logdbskewncopu <- function(u, v, cpar) logdbskewncop(1-u, v, cpar)

## Conditional distribution

#' @rdname bskewncop_supp
#' @export
pcondbskewncopu <- function(v, u, cpar) pcondbskewncop(v, 1-u, cpar)
#' @rdname bskewncop_supp
#' @export
qcondbskewncopu <- function(p, u, cpar) qcondbskewncop(p, 1-u, cpar)

## ----- V-Flipped / U-axis reflection / Vertical reflection -----

## Distribution

#' @rdname bskewncop_supp
#' @export
pbskewncopv <- function(u, v, cpar) u - pbskewncop(u, 1-v, cpar)

## Density

#' @rdname bskewncop_supp
#' @export
dbskewncopv <- function(u, v, cpar) dbskewncop(u, 1-v, cpar)
#' @rdname bskewncop_supp
#' @export
logdbskewncopv <- function(u, v, cpar) logdbskewncop(u, 1-v, cpar)

## Conditional distribution

#' @rdname bskewncop_supp
#' @export
pcondbskewncopv <- function(v, u, cpar) 1 - pcondbskewncop(1-v, u, cpar)
#' @rdname bskewncop_supp
#' @export
qcondbskewncopv <- function(p, u, cpar) 1 - qcondbskewncop(1-p, u, cpar)

