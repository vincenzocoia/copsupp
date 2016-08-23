#' Construction Function for the IG Copula Family
#'
#' \code{cnstr_H} is the function itself, and \code{cnstr_Hinv} is
#' its inverse; \code{cnstr_D1H} is the first-argument
#' derivative.
#'
#' @param t Vector of values >=1 to evaluate the function at.
#' @param theta Value of second argument of H_{k}, >0. This is allowed
#' to be a vector, except in \code{cnstr_Hinv()}.
#' @param k Single numeric >1 corresponding to the \code{k} parameter
#' of the function.
#' @rdname cnstr_H
#' @export
cnstr_H <- function(t, theta, k) {
    cnstr_Psi(1/(theta * log(t)), k) / t
}

#' @rdname cnstr_H
#' @export
cnstr_D1H <- function(t, theta, k) {
    fun <- function(t) {
        logt <- log(t)
        arg <- 1/theta/logt
        coeff <- 1/theta/logt^2
        -t^(-2) * (cnstr_Psi(arg,k) + coeff * cnstr_DPsi(arg,k))
    }
    ## Deal with t=1 separately -- its limit depends on k.
    if (k > 2) {
        replf <- -1
    } else if (k == 2) {
        replf <- -(1 + theta/2)
    } else {
        replf <- -Inf
    }
    eval_lims(fun, t, replx=1, replf=replf)
}


#' @param w Vector of values in [0,1] to evaluate the inverse function at.
#' @param silent Logical; should the message output by \code{pcinterpolate()}
#' be silenced?
#' @rdname cnstr_H
#' @import CopulaModel
#' @export
cnstr_Hinv <- function(w, theta, k, ngrid=1000, silent=TRUE) {
    ## Main Idea: -----
    ## Get values to evaluate cnstr_H at, by using an approximation function
    ##  that can be inverted. We'll invert that function on a grid in (0,1)
    ##  to get a grid on the domain of H.
    ## It seems to be easier to find an approximation function by transforming
    ##  the domain from (1,oo) to (-oo,oo) by taking a log transform twice,
    ##  so we'll do that.
    ## ----------------
    fun <- function(clean_w) {
        ## Get grid on domain of H:
        wmin <- min(0.001, min(clean_w)/2)
        wmax <- max(0.999, (1+max(clean_w))/2)
        wgrid <- seq(wmin, wmax, length.out=ngrid)
        tgrid <- -log(-log(1-wgrid)) - log(1 + theta/k)
        ## Evaluate H at that grid:
        fn <- cnstr_H(exp(exp(tgrid)), theta, k)
        ## The evaluated cnstr_H become the domain values of cnstr_Hinv, and
        ##  the tgrid values become the evaluated cnstr_Hinv function.
        ## Use the pcinterpolate function to evaluate at w.
        der <- pcderiv(fn, tgrid)
        ## Evaluate at the requested NA-free w's
        if (silent) {
            ## The only way I can find to silent a cat() call is to use sink,
            ##  but send it to a non-existing file. Got answer from Stack
            ##  Overflow by cbielow:
            ## http://stackoverflow.com/questions/6177629/how-to-silence-the-output-from-this-r-package
            f <- file()
            sink(file=f)
            res <- pcinterpolate(fn, tgrid, der, clean_w)[, 1]
            sink()
            close(f)
        } else {
            res <- pcinterpolate(fn, tgrid, der, clean_w)[, 1]
        }
        ## Put back onto (1,oo) domain:
        return(exp(exp(res)))
    }
    eval_lims(fun, w, replx=c(0,1), replf=c(Inf, 1))
}
