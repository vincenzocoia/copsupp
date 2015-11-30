#' Multivariate Integration
#'
#' Integrate a multivariate function (R^p -> R).
#'
#' @param f Function that takes a vector and returns a numeric.
#' @param lower, upper Vector of the limits of integration for
#' each variable entered into \code{f}. Can be infinite.
#' @param stop.on.error Logical (an argument of \code{\link{integrate}}) --
#' should function be stopped if an "error" occurs?
#' @param ... Other arguments to pass into \code{\link{integrate}}.
#' @details This function recursively integrates the arguments in \code{f}
#' using the \code{\link{integrate}} function.
#'
#' It's reasonably fast
#' when \code{f} takes two variables, and slow when \code{f} takes 3 variables.
#' For 4 or more, perhaps look into alternative approaches.
#'
#' If you're support is finite, you're better off using
#' \code{\link{cubature::adaptIntegrate}}.
#' @return A numeric value of the integral. The other information
#' outputted by the \code{\link{integrate}} function are ignored.
#'
#' If the length of the limit vectors is zero, it will be assumed that
#' \code{f} has no arguments, so \code{f()} will be returned.
#' returned.
#' @note Only a rectangular support is allowed.
#'
#' This function could probably be improved if some of the support is finite,
#' in which case somehow \code{\link{cubature::adaptIntegrate}} could be
#' leveraged somehow.
#'
#' The \code{stop.on.error} argument is defaulted to \code{FALSE}
#' instead of \code{\link{integrate}}'s \code{TRUE} because sometimes
#' it'll think the integral is divergent when it's really not.
#' @examples
#' pdf <- function(x) prod(dnorm(x))
#' integrate.mv(pdf, -Inf, Inf)
#' integrate.mv(pdf, c(-Inf, -Inf), c(Inf, Inf))
#' @export
integrate.mv <- function(f, lower, upper, stop.on.error = FALSE, ...) {
    p <- length(lower)
    if (p != length(upper)) stop("Lower and upper endpoints of unequal length.")
    if (p == 0) return(f())
    if (p == 1) {
        integrand <- Vectorize(f)
        return(integrate(integrand, lower, upper, stop.on.error = stop.on.error, ...)$value)
    } else {
        nextf <- function(firstargs) {
            integrand_ <- function(xp) f(c(firstargs, xp))
            integrand <- Vectorize(integrand_)
            integrate(integrand, lower[p], upper[p], stop.on.error = stop.on.error, ...)$value
        }
        return(integrate.mv(nextf, lower[-p], upper[-p], stop.on.error = stop.on.error, ...))
    }
}
