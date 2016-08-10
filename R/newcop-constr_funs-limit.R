#' Construction functions for the new copula
#' @examples
#' ell(0:10, k=3)
#' ell(0:10, k=1)
#' ell(0:10, k=99999)
#' @rdname const_fun
#' @export
cnstr_ell <- function(t, k) {
    ## t's that are zero should evaluate to 1.
    whichzeroes <- t==0
    res <- rep(NA, length(t))
    res[whichzeroes] <- 1
    tt <- t[!whichzeroes]
    ## non-zero t's should evaluate according to the formula.
    res[!whichzeroes] <- 1 - pgamma(tt, k-1) + (k-1)/tt*pgamma(tt, k)
    res
}

#' @rdname const_fun
#' @export
cnstr_ellp <- function(k, t) {
    -(k-1)/t^2 * pgamma(t, k)
}

#' @rdname const_fun
#' @export
cnstr_ellp2 <- function(k, t) {
    (k-1) * (2*(k-1)/t^3 * pgamma(t,k) - t^(k-3)*exp(-t)/gamma(k-1))
}
