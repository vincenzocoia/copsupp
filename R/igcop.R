#' IG Copula Family Functions
#'
#' Functions related to the IG copula family, denoted  by \code{'igcop'}.
#'
#' @param u,v Vectors of values in [0,1] representing values of the first
#' and second copula variables.
#' @param cpar Vector of length 2 corresponding to the copula
#' parameters \code{theta>0} and \code{k>1}, respectively.
#' @note Inputting two vectors greater than length 1 is allowed, if they're
#' the same length.
#' Also, \code{qcondigcop21} and \code{pcondigcop21} are the same as
#' \code{qcondigcop} and \code{pcondigcop} -- their the distributions of
#' variable 2 given 1.
#' @return Numeric vector of length equal to the length of the input vector(s).
#' @rdname igcop
#' @export
pcondigcop <- function(v, u, cpar) {
    theta <- cpar[1]
    k <- cpar[2]
    Hkinv <- cnstr_Hinv(1-v, theta, k)
    1 - pgamma(theta * (1-u) * log(Hkinv), k-1, lower.tail=FALSE) / Hkinv
}



#' @param tau Vector of quantile levels in [0,1] to evaluate a quantile function
#' at.
#' @note Of all the methods of computing qcondigcop() that I tried,
#' this is the only one that appears to work (validated
#' visually with the old qcondnew()
#' function for theta=3 and 30). This version uses the inverse of the
#' helper function, \code{\link{igcop_helper_inv}}, but deliberately
#' inefficiently by evaluating each helper function at all the taus, just
#' because the inverse seems to work better that way.
#' @rdname igcop
#' @export
qcondigcop <- function(tau, u, cpar) {
    ## Get "parameters" useful for this function.
    theta <- cpar[1]
    k <- cpar[2]
    param <- (1 - u) * theta
    n_u <- length(u)
    n_tau <- length(tau)
    ## Was only one tau input? If so, repeat it to match the length of u.
    if (n_tau == 1) tau <- rep(tau, n_u)
    ## If u is of length 1, then we can use vectorized property of
    ##  igcop_helper_inv().
    if (n_u == 1) {
        fun <- function(tau) {
            arg <- exp(exp(igcop_helper_inv(tau, param, k)))
            1 - cnstr_H(arg, theta, k)
        }
        eval_lims(fun, tau, replx=c(0,1), replf=c(0,1))
    } else {
        ## There's a cdf for each u that needs inverting. Loop through and
        ##  invert them one-by-one.
        mapply(function(tau_, param_) {
            fun <- function(tau_) {
                arg <- exp(exp(igcop_helper_inv(tau_, param_, k)))
                1 - cnstr_H(arg, theta, k)
            }
            eval_lims(fun, tau_, replx=c(0,1), replf=c(0,1))
        }, tau, param)
    }
}

#' @rdname igcop
#' @export
digcop <- function(u, v, cpar) {
    theta <- cpar[1]
    k <- cpar[2]
    negu <- 1-u
    t <- cnstr_Hinv(1-v, theta, k)
    x <- theta * negu * log(t)
    -(dgamma(x, k-1) * theta + pgamma(x, k-1, lower.tail=FALSE)) / t^2 /
        cnstr_D1H(t, theta, k)
}
