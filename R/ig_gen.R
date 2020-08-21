#' Generating Function for the IG Copula Family
#'
#' \code{ig_gen} is the function itself, and \code{ig_geninv} is
#' its inverse; \code{ig_D1gen} is the first-argument
#' derivative.
#'
#' @param t Vector of values >=1 to evaluate the function at.
#' @param eta Value of second argument of H_{k}, >0. This is allowed
#' to be a vector, except in \code{ig_geninv()}.
#' @param k Single numeric >1 corresponding to the \code{k} parameter
#' of the function.
#' @rdname ig_gen
#' @export
ig_gen <- function(t, eta, k) {
    igl_gen(1/(eta * log(t)), k) / t
}

#' @rdname ig_gen
#' @export
ig_D1gen <- function(t, eta, k) {
    fun <- function(t) {
        logt <- log(t)
        arg <- 1/eta/logt
        coeff <- 1/eta/logt^2
        -t^(-2) * (igl_gen(arg,k) + coeff * igl_Dgen(arg,k))
    }
    ## Deal with t=1 separately -- its limit depends on k.
    if (k > 2) {
        replf <- -1
    } else if (k == 2) {
        replf <- -(1 + eta/2)
    } else {
        replf <- -Inf
    }
    eval_lims(fun, t, replx=1, replf=replf)
}


#' @param w Vector of values in [0,1] to evaluate the inverse function at.
#' @rdname ig_gen
#' @import CopulaModel
#' @export
ig_geninv <- function(p, eta, k, mxiter=40, eps=1.e-12, bd=5) {
    ## Algorithm:
    fun <- function(p, eta) {
        ## Get starting values
        xp1 <- 1/p
        xp2 <- exp(1/eta/igl_geninv(p, k))
        xpm <- pmin(xp1,xp2)
        tt <- pmax(xpm - eps, 1 + (xpm-1)/2) # xpm-eps might overshoot left of 1.
        tt <- pmax(1 + eps, tt)  # tt might be 1.
        iter <- 0
        diff <- 1
        ## Begin Newton-Raphson algorithm
        while(iter<mxiter & max(abs(diff))>eps){
            ## Helpful quantities
            logt <- log(tt)
            etalog <- eta * logt
            etaloginv <- 1 / etalog
            ## Evaluate functions
            g <- tt*p - igl_gen(etaloginv, k)
            gp <- p + etaloginv/tt/logt * igl_Dgen(etaloginv, k)
            diff <- g/gp
            flag <- diff > tt - 1
            diff[flag] <- (tt[flag] - 1) / 2
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
