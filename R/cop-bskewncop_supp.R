## ----- Permutation-reflected -----

## Distribution

#' Skew normal copula: supplement
#'
#' @rdname bskewncop_supp
#' @export
pbskewncopp <- function(u, v, cpar) pbskewncop(v, u, cpar)

## Density

#' @rdname bskewncop_supp
#' @export
dbskewncopp <- function(u, v, cpar) dbskewncop(v, u, cpar)
#' @rdname bskewncop_supp
#' @export
logdbskewncopv <- function(u, v, cpar) logdbskewncop(v, u, cpar)

## Conditional distribution

#' @rdname bskewncop_supp
#' @export
pcondbskewncopp21 <- function(v, u, cpar) pcondbskewncop12(v, u, cpar)
#' @rdname bskewncop_supp
#' @export
qcondbskewncopp21 <- function(p, u, cpar) qcondbskewncopp12(p, u, cpar)
#' @rdname bskewncop_supp
#' @export
pcondbskewncopp <- function(v, u, cpar) pcondbskewncop12(v, u, cpar)
#' @rdname bskewncop_supp
#' @export
qcondbskewncopp <- function(p, u, cpar) qcondbskewncopp12(p, u, cpar)
#' @rdname bskewncop_supp
#' @export
pcondbskewncopp12 <- function(u, v, cpar) pcondbskewncop21(u, v, cpar)
#' @rdname bskewncop_supp
#' @export
qcondbskewncopp12 <- function(p, v, cpar) qcondbskewncopp21(p, v, cpar)


## ---- Regular ----

#' @rdname bskewncop_supp
#' @export
pcondbskewncop <- function(v, u, cpar) pcondbskewncop21(v, u, cpar)
#' @rdname bskewncop_supp
#' @export
qcondbskewncop <- function(p, u, cpar) qcondbskewncop21(p, u, cpar)
