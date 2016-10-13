#' Functions for IG copula 2|1 computations
#'
#' @param t Vector of values to evaluate the assistant function at, >=1.
#' @param p Vector of values to evaluate the inverse assistant function at,
#' in [0,1].
#' @param k Single numeric, greater than 1.
#' @param eta eta parameter. Vectorized to match \code{t} and \code{p}.
#' @param mxiter Maximum number of iterations to run the Newton-Raphson
#' algorithm for.
#' @param eps The Newton-Raphson
#' algorithm is stopped if all step sizes are below this value.
#' @param bd Largest acceptable step size in the Newton-Raphson algorithm;
#' steps larger than this will be reduced.
#' @rdname igcond
#' @export
igcond <- function(t, k, eta) {
    fun <- function(t, eta) pgamma(eta*log(t), k-1, lower.tail=FALSE) / t
    eval_lims(fun, t, replx=Inf, replf=0, eta=eta)
}

#' @rdname igcond
#' @export
igcondinv <- function(p, k, eta, mxiter=40, eps=1.e-6, bd=5) {
    ## Algorithm:
    fun <- function(p, eta) {
        ## Get starting values
        xp1 <- 1/p
        xp2 <- exp(qgamma(1-p,k-1)/eta)
        xpm <- pmin(xp1,xp2)
        tt <- pmax(xpm - eps, 1 + (xpm-1)/2) # xpm-eps might overshoot left of 1.
        iter <- 0
        diff <- 1
        ## Begin Newton-Raphson algorithm
        while(iter<mxiter & max(abs(diff))>eps){
            ## Helpful quantities
            etalog <- eta * log(tt)
            ## Evaluate functions
            g <- tt * p - pgamma(etalog, k-1, lower.tail=FALSE)
            ## When eta=0, derivative is NaN when 1<k<2, when should just be p.
            gpfun <- function(eta) p + dgamma(etalog, k-1) * eta / tt
            gp <- eval_lims(gpfun, eta, replx=0, replf=p)
            diff <- g/gp
            tt <- tt-diff
            while(max(abs(diff))>bd | any(tt<=1))
            { diff <- diff/2; tt <- tt+diff }
            iter <- iter+1
            # cat(paste0("-----", iter, "-----\n"))
            # cat(diff, "\n")
            # cat(tt, "\n")
        }
        return(tt)
    }
    eval_lims(fun, p, replx=c(0, 1), replf=c(Inf, 1), eta=eta)
}


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
    if (theta == Inf) return(pcondiglcop(v, u, k))
    Hkinv <- cnstr_Hinv(1-v, theta, k)
    1 - igcond(Hkinv, k, theta * (1-u))
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
    theta <- cpar[1]
    k <- cpar[2]
    if (theta == Inf) return(qcondiglcop(tau, u, k))
    inv <- igcondinv(1-tau, k, theta*(1-u))
    1 - cnstr_H(inv, theta, k)
}

#' @rdname igcop
#' @export
digcop <- function(u, v, cpar) {
    theta <- cpar[1]
    k <- cpar[2]
    if (theta == Inf) return(diglcop(u, v, k))
    negu <- 1-u
    t <- cnstr_Hinv(1-v, theta, k)
    x <- theta * negu * log(t)
    -(dgamma(x, k-1) * theta + pgamma(x, k-1, lower.tail=FALSE)) / t^2 /
        cnstr_D1H(t, theta, k)
}

#' @rdname igcop
#' @export
logdigcop <- function(u, v, cpar) log(digcop(u, v, cpar))

#' @rdname igcop
#' @export
pigcop <- function(u, v, cpar) {
    theta <- cpar[1]
    k <- cpar[2]
    Hinv <- cnstr_Hinv(1-v, theta, k)
    u + v - 1 + (1-u) * cnstr_H(Hinv, theta * (1-u))
}
