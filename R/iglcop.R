#' IGL Copula Family Functions
#'
#' Functions related to the IGL copula family, denoted  by \code{'iglcop'}.
#'
#' @param u,v Vectors of values in [0,1] representing values of the first
#' and second copula variables.
#' @param tau Vector of quantile levels in [0,1] to evaluate a quantile function
#' at.
#' @param cpar Single numeric >1; corresponds to parameter \code{k} in the
#' IGL copula family.
#' @note Inputting two vectors greater than length 1 is allowed, if they're
#' the same length.
#' Also, \code{qcondiglcop21} and \code{pcondiglcop21} are the same as
#' \code{qcondiglcop} and \code{pcondiglcop} -- their the distributions of
#' variable 2 given 1.
#' @return Numeric vector of length equal to the length of the input vector(s).
#' @rdname iglcop
#' @export
qcondiglcop <- function(tau, u, cpar) {
    1 - igl_gen((1 - u)/qgamma(tau, cpar-1), cpar)
}

#' @rdname iglcop
#' @export
pcondiglcop <- function(v, u, cpar) {
    pgamma((1 - u) / igl_geninv(1 - v, cpar), cpar - 1)
}

#' @rdname iglcop
#' @export
qcondiglcop21 <- qcondiglcop

#' @rdname iglcop
#' @export
pcondiglcop21 <- pcondiglcop

#' @rdname iglcop
#' @export
pcondiglcop12 <- function(u, v, cpar) {
    pkinv <- igl_geninv(1 - v, cpar)
    1 - pgamma((1 - u) / pkinv, cpar) / pgamma(1 / pkinv, cpar)
}

#' @rdname iglcop
#' @export
qcondiglcop12 <- function(tau, v, cpar) {
    pkinv <- igl_geninv(1 - v, cpar)
    1 - pkinv * qgamma((1 - tau) * pgamma(1/pkinv, cpar), cpar)
}

#' @rdname iglcop
#' @export
diglcop <- function(u, v, cpar) {
    pkinv <- igl_geninv(1 - v, cpar)
    (1-u)^(cpar-1) / pkinv^cpar * exp(-(1-u)/pkinv) / (gamma(cpar) * pgamma(1/pkinv, cpar))
}

#' @rdname iglcop
#' @export
piglcop <- function(u, v, cpar) {
    pkinv <- igl_geninv(1-v, cpar)
    u + v - 1 + (1-u) * igl_gen(pkinv / (1-u), cpar)
}

#' @rdname igcop
#' @export
riglcop <- function(n, cpar) {
    u <- runif(n)
    tau <- runif(n)
    v <- qcondiglcop(tau, u, cpar)
    matrix(c(u, v), ncol=2)
}
