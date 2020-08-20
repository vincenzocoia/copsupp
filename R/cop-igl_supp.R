

## ----- U-Flipped / V-axis reflection / Horizontal reflection -----

## Distribution

#' IG copula: supplemental functions
#'
#' Appended "u" means a horizontal flip (i.e. switching u to 1-u).
#' Appended "v" means a vertical flip (i.e. switching v to 1-v).
#' Appended "r" means reflection/survival copula
#'
#' @rdname iglcop_supp
#' @export
piglcopu <- function(u, v, cpar) v - piglcop(1-u, v, cpar)

## Density

#' @rdname iglcop_supp
#' @export
diglcopu <- function(u, v, cpar) diglcop(1-u, v, cpar)
#' @rdname iglcop_supp
#' @export
logdiglcopu <- function(u, v, cpar) logdiglcop(1-u, v, cpar)

## Random number generator

#' @rdname iglcop_supp
#' @export
riglcopu <- function(n, cpar) {
    res <- riglcop(n, cpar)
    res[, 1] <- 1 - res[, 1]
    res
}

## Conditional distribution

#' @rdname iglcop_supp
#' @export
pcondiglcopu <- function(v, u, cpar) pcondiglcop(v, 1-u, cpar)
#' @rdname iglcop_supp
#' @export
qcondiglcopu <- function(p, u, cpar) qcondiglcop(p, 1-u, cpar)

## ----- V-Flipped / U-axis reflection / Vertical reflection -----

## Distribution

#' @rdname iglcop_supp
#' @export
piglcopv <- function(u, v, cpar) u - piglcop(u, 1-v, cpar)

## Density

#' @rdname iglcop_supp
#' @export
diglcopv <- function(u, v, cpar) diglcop(u, 1-v, cpar)
#' @rdname iglcop_supp
#' @export
logdiglcopv <- function(u, v, cpar) logdiglcop(u, 1-v, cpar)

## Conditional distribution

#' @rdname iglcop_supp
#' @export
pcondiglcopv <- function(v, u, cpar) 1 - pcondiglcop(1-v, u, cpar)
#' @rdname iglcop_supp
#' @export
qcondiglcopv <- function(p, u, cpar) 1 - qcondiglcop(1-p, u, cpar)


## ----- Reflection ------


#' @rdname iglcop_supp
#' @export
pcondiglcopr <- function(v, u, cpar) 1 - pcondiglcop(1 - v, 1 - u, cpar)

#' @rdname iglcop_supp
#' @export
qcondiglcopr <- function(p, u, cpar) 1 - qcondiglcop(1 - p, 1 - u, cpar)

#' @rdname iglcop_supp
#' @export
pcondiglcopr12 <- function(u, v, cpar) 1 - pcondiglcop12(1 - u, 1 - v, cpar)

#' @rdname iglcop_supp
#' @export
qcondiglcopr12 <- function(p, v, cpar) 1 - qcondiglcop12(1 - p, 1 - v, cpar)

#' @rdname iglcop_supp
#' @export
diglcopr <- function(u, v, cpar) diglcop(1 - u, 1 - v, cpar)

#' @rdname iglcop_supp
#' @export
piglcopr <- function(u, v, cpar) u + v - 1 + piglcop(1 - u, 1 - v, cpar)

#' @rdname iglcop_supp
#' @export
riglcopr <- function(n, cpar) 1 - riglcop(n, cpar)
