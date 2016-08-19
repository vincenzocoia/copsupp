#' IG Copula Family Functions
#'
#' Functions related to the IG copula family, denoted  by \code{'igcop'}.
#'
#' @param u,v Vectors of values in [0,1] representing values of the first
#' and second copula variables.
#' @param tau Vector of quantile levels in [0,1] to evaluate a quantile function
#' at.
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
        arg <- mapply(function(tau_, param_) {
            exp(exp(igcop_helper_inv(tau_, param_, k)))
        }, clean_tau, clean_param)
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


# qcondigcop <- function(tau, u, cpar, mxiter=20,eps=1.e-6,bd=5){
#     ## Make tau and u of the same length
#     n_tau <- length(tau)
#     n_u <- length(u)
#     n <- max(n_tau, n_u)
#     if (n_tau == 1) tau <- rep(tau, n)
#     if (n_u == 1) u <- rep(u, n)
#     # ## Work with non-NA, non-1, non-0 values.
#     # NAs_tau <- is.na(tau)
#     # NAs_u <- is.na(u)
#     # NAs <-
#     # ones <- (tau == 1)  # T/F. Has NA's too.
#     # whichones <- which(ones)
#     # zeroes <- (tau == 0) # T/F. Has NA's too.
#     # whichzeroes <- which(zeroes)
#     # clean_tau <- na.omit(tau[!(ones | zeroes)])
#     clean_tau <- tau
#     ## Go ahead with the algorithm
#     theta <- cpar[1]
#     k <- cpar[2]
#     if (length(clean_tau) > 0){
#         ## Use halfway between independence and comonotonicity.
#         tt <- (clean_tau + u)/2
#         iter <- 0
#         diff <- 1
#         ## Begin Newton-Raphson algorithm
#         while(iter<mxiter & max(abs(diff))>eps){
#             ## Helpful quantities
#             Hkinv <- cnstr_Hinv(1-tt, theta, k)
#             arg <- theta * (1-u) * log(Hkinv)
#             der <- cnstr_D1H(Hkinv, theta, k)
#             dgam <- dgamma(arg, k-1) * theta * (1-u)
#             ## Evaluate functions
#             g <- dgam/der/Hkinv + (1-tau)/der
#             gp <- 1/der * (dgam/Hkinv + 1 - tau)
#             diff <- g/gp
#             tt <- tt-diff
#             while(max(abs(diff))>bd | any(tt<=0))
#             { diff <- diff/2; tt <- tt+diff }
#             iter <- iter+1
#             #cat(iter,diff,tt,"\n")
#         }
#     } else {
#         tt <- numeric(0)
#     }
#
#     # ## Set up vector to be returned (start off with NA's)
#     # res <- rep(NA, max(length(tau), length(u)))
#     # special_indices <- c(which(NAs), whichones, whichzeroes)
#     # if (length(special_indices > 0)){
#     #     res[-special_indices] <- tt
#     #     res[whichones] <- 1
#     #     res[whichzeroes] <- 0
#     # } else {
#     #     res <- tt
#     # }
#     # res
#     return(tt)
# }
