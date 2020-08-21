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
    Hkinv <- ig_geninv(1-v, theta, k)
    1 - igcond(Hkinv, k, theta * (1-u))
}

#' @param tau Vector of quantile levels in [0,1] to evaluate a quantile function
#' at.
#' @rdname igcop
#' @export
qcondigcop <- function(tau, u, cpar) {
    theta <- cpar[1]
    k <- cpar[2]
    if (theta == Inf) return(qcondiglcop(tau, u, k))
    inv <- igcondinv(1-tau, k, theta*(1-u))
    1 - ig_gen(inv, theta, k)
}

#' @rdname igcop
#' @export
qcondigcop21 <- qcondigcop

#' @rdname igcop
#' @export
pcondigcop21 <- pcondigcop

#' @rdname igcop
#' @export
pcondigcop12 <- function(u, v, cpar) {
    theta <- cpar[1]
    k <- cpar[2]
    if (theta == Inf) return(pcondiglcop12(u, v, k))
    Hkinv <- ig_geninv(1-v, theta, k)
    1 - (1-u) * ig_D1gen(Hkinv, theta*(1-u), k) / ig_D1gen(Hkinv, theta, k)
}

#' @rdname igcop
#' @export
digcop <- function(u, v, cpar) {
    theta <- cpar[1]
    k <- cpar[2]
    if (theta == Inf) return(diglcop(u, v, k))
    # negu <- 1-u
    # t <- ig_geninv(1-v, theta, k)
    # x <- theta * negu * log(t)
    # -(dgamma(x, k-1) * theta + pgamma(x, k-1, lower.tail=FALSE)) / t^2 /
    #     ig_D1gen(t, theta, k)
    u <- 1 - u
    v <- 1 - v
    tv <- ig_geninv(v, theta, k)
    ltv <- log(tv)
    # num1 <- (theta * u) ^ (k - 1) * ltv ^ (k - 2) * tv ^ (-theta * u) / gamma(k - 1)
    num1 <- theta * u * stats::dgamma(theta * u * ltv, k - 1)
    num2 <- stats::pgamma(theta * u * ltv, k - 1, lower.tail = FALSE)
    den1 <- (ltv + 1) * (k - 1) / (theta * ltv ^ 2) * stats::pgamma(theta * ltv, k)
    den2 <- stats::pgamma(theta * ltv, k - 1, lower.tail = FALSE)
    (num1 + num2) / (den1 + den2)
}

#' @rdname igcop
#' @export
logdigcop <- function(u, v, cpar) log(digcop(u, v, cpar))

#' @rdname igcop
#' @export
pigcop <- function(u, v, cpar) {
    theta <- cpar[1]
    k <- cpar[2]
    Hinv <- ig_geninv(1-v, theta, k)
    u + v - 1 + (1-u) * ig_gen(Hinv, theta * (1-u), k)
}

#' @rdname igcop
#' @export
rigcop <- function(n, cpar) {
    u <- runif(n)
    tau <- runif(n)
    v <- qcondigcop(tau, u, cpar)
    matrix(c(u, v), ncol=2)
}
