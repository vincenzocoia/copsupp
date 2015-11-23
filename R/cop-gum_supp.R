## ----- U-Flipped / V-axis reflection / Horizontal reflection -----

## Distribution

#' Gumbel copula: supplemental functions
#'
#' Appended "u" means a horizontal flip (i.e. switching u to 1-u).
#' Appended "v" means a vertical flip (i.e. switching v to 1-v).
#' Appended "r" means reflection/survival copula
#'
#' @rdname gum_supp
#' @export
pgumu <- function(u, v, cpar) v - pgum(1-u, v, cpar)

## Density

#' @rdname gum_supp
#' @export
dgumu <- function(u, v, cpar) dgum(1-u, v, cpar)
#' @rdname gum_supp
#' @export
logdgumu <- function(u, v, cpar) logdgum(1-u, v, cpar)

## Random number generator

#' @rdname gum_supp
#' @export
rgumu <- function(n, cpar) {
    res <- rgum(n, cpar)
    res[, 1] <- 1 - res[, 1]
    res
}

## Conditional distribution

#' @rdname gum_supp
#' @export
pcondgumu <- function(v, u, cpar) pcondgum(v, 1-u, cpar)
#' @rdname gum_supp
#' @export
qcondgumu <- function(p, u, cpar) qcondgum(p, 1-u, cpar)

## ----- V-Flipped / U-axis reflection / Vertical reflection -----

## Distribution

#' @rdname gum_supp
#' @export
pgumv <- function(u, v, cpar) u - pgum(u, 1-v, cpar)

## Density

#' @rdname gum_supp
#' @export
dgumv <- function(u, v, cpar) dgum(u, 1-v, cpar)
#' @rdname gum_supp
#' @export
logdgumv <- function(u, v, cpar) logdgum(u, 1-v, cpar)

## Random number generator

#' @rdname gum_supp
#' @export
rgumv <- function(n, cpar) {
    res <- rgum(n, cpar)
    res[, 2] <- 1 - res[, 2]
    res
}

## Conditional distribution

#' @rdname gum_supp
#' @export
pcondgumv <- function(v, u, cpar) 1 - pcondgum(1-v, u, cpar)
#' @rdname gum_supp
#' @export
qcondgumv <- function(p, u, cpar) 1 - qcondgum(1-p, u, cpar)

