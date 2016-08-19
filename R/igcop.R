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
    ## Which taus are 0 and 1? Treat them specially. u is allowed to be 1,
    ##  because igcop_helper_inv can handle theta=0. u=0 results in evaluating
    ##  its theta parameter at (1-u)*theta=theta, so no problem there.
    zeroes <- tau == 0  # Keeps NAs
    ones <- tau == 1  # Keeps NAs
    whichzeroes <- which(zeroes)  # NAs not included.
    whichones <- which(ones)  # NAs not included.
    ## Which are NA? This goes for both tau and u!
    NAs <- is.na(tau) | is.na(u)
    which_NA <- which(NAs)
    ## Only consider non-NA values of tau in (0,1). Even if this subsets
    ##  to an empty vector, the following procedure should still work.
    rmv <- union(which_NA, c(whichzeroes, whichones))
    if (length(rmv) > 0) {
        clean_tau <- tau[-rmv]  # Note tau[-integer(0)] = numeric(0)
    } else {
        clean_tau <- tau
    }
    ## If u is of length 1, then we can use vectorized property of
    ##  igcop_helper_inv().
    if (n_u == 1) {
        arg <- exp(exp(igcop_helper_inv(clean_tau, param, k)))
        res1 <- 1 - cnstr_H(arg, theta, k)
    } else {
        ## u (and param) is not of length 1. Suppose it's of length tau.
        ## Now subset param to match those of clean_tau.
        if (length(rmv) > 0) {
            clean_param <- param[-rmv]  # Note tau[-integer(0)] = numeric(0)
        } else {
            clean_param <- param
        }
        ## Get inverse helper function for each param value.
        arg <- sapply(1:length(clean_tau), function(i){
            exp(exp(igcop_helper_inv(clean_tau, clean_param[i], k)[i]))
        })
        ## And the "clean" quantile values:
        res1 <- 1 - cnstr_H(arg, theta, k)
    }
    ## Lastly, put NAs where they belong, and replace tau=0 and tau=1
    ##  with quantiles of 0 and 1.
    if (length(rmv) == 0) {
        return(res1)
    } else {
        res <- rep(NA, length(tau))
        res[whichzeroes] <- 0
        res[whichones] <- 1
        res[-rmv] <- res1
        return(res)
    }
}
