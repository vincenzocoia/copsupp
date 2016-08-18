#' Construction function for the IGL copula family
#'
#' \code{cnstr_Psi} is the function itself, and \code{cnstr_Psiinv} is
#' its inverse; \code{cnstr_DPsi} is the
#' derivative; and \code{cnstr_D2Psi} is the second derivative.
#'
#' @param t Vector of values to evaluate the function at, >=0.
#' @param w Vector of values to evaluate the inverse function at, between
#' 0 and 1 (inclusive)
#' @param k Single numeric >1. Parameter of the function.
#' @examples
#' ## Some examples of evaluating the functions.
#' arg <- c(0, 0.5, 3, Inf, NA)
#' cnstr_Psi(arg, k=2)
#' cnstr_DPsi(arg, k=1.2)
#' cnstr_DPsi(arg, k=2)
#' cnstr_DPsi(arg, k=3)
#' cnstr_Psiinv(c(0, 0.5, 1), k=1.5)
#'
#' ## Visual
#' foo <- function(u) cnstr_Psiinv(u, k=1.5)
#' curve(foo)
#' @rdname cnstr_Psi
#' @export
cnstr_Psi <- function(t, k) {
    ## Allow for the possibility of t=Inf (this needs special attention
    ##   because Inf*pgamma(0,k) = NaN)
    whichInf <- which(t == Inf)  # NA's won't appear in the vector, because of which().
    whichelse <- setdiff(seq_len(length(t)), whichInf)
    ## Note: whichelse is needed, because t[-integer(0)] does not equal t -- in
    ##  fact, it's empty! Access non-infinity values by t[whichelse].
    ## Set up the result vector
    res <- rep(NA, length(t))
    ## Compute the function at non-infinite values:
    tt <- t[whichelse]
    tinv <- 1/tt
    res[whichelse] <- 1 - pgamma(tinv, k-1) + (k-1)*tt*pgamma(tinv, k)
    ## At infinity, the function is 1. Put that in the output
    res[whichInf] <- 1
    return(res)
}

#' @rdname cnstr_Psi
#' @export
cnstr_DPsi <- function(t, k) {
    (k-1) * pgamma(1/t, k)
}

#' @rdname cnstr_Psi
#' @export
cnstr_D2Psi <- function(t, k) {
    -t^(-k-1) * exp(-1/t) / gamma(k-1)
}

#' @rdname cnstr_Psi
#' @export
cnstr_Psiinv <- function(w, k, mxiter=20,eps=1.e-6,bd=5){
    ## Work with non-NA, non-1, non-0 values.
    NAs <- is.na(w)
    ones <- (w == 1)  # T/F. Has NA's too.
    whichones <- which(ones)
    zeroes <- (w == 0) # T/F. Has NA's too.
    whichzeroes <- which(zeroes)
    clean_w <- na.omit(w[!(ones | zeroes)])
    ## Compute gamma(k-1) and gamma(k)
    gkm1 <- gamma(k-1)
    gk <- (k-1) * gkm1
    ## Go ahead with the algorithm
    if (length(clean_w) > 0){
        v <- clean_w
        ## Empirically, it looks like this tt is a good start:
        tt <- (1-v)^(-1/(k-1)) - 1
        iter <- 0
        diff <- 1
        ## Begin Newton-Raphson algorithm
        while(iter<mxiter & max(abs(diff))>eps){
            ## Helpful quantities
            igam1 <- gkm1 * pgamma(1/tt, k-1, lower.tail=FALSE)
            igam0 <- gk * pgamma(1/tt, k)
            ## Evaluate functions
            g <- tt * igam0 + igam1 - v * gkm1
            gp <- igam0
            diff <- g/gp
            tt <- tt-diff
            while(max(abs(diff))>bd | any(tt<=0))
            { diff <- diff/2; tt <- tt+diff }
            iter <- iter+1
            #cat(iter,diff,tt,"\n")
        }
    } else {
        tt <- numeric(0)
    }

    ## Set up vector to be returned (start off with NA's)
    res <- rep(NA, length(w))
    special_indices <- c(which(NAs), whichones, whichzeroes)
    if (length(special_indices > 0)){
        res[-special_indices] <- tt
        res[whichones] <- Inf
        res[whichzeroes] <- 0
    } else {
        res <- tt
    }
    res
}
