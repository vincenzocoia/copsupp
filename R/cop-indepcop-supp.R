#' Conditional Independence Copula
#'
#' These are just the identity functions with respect to the first argument.
#' @rdname indepcopcond
#' @export
pcondindepcop <- function(v,u,cpar=0){
    nu <- length(u)
    nv <- length(v)
    if (nv == 1) {
        return(rep(v, nu))
    } else {
        if (nu != nv & nu != 1) {
            stop("u and v must be the same length, or one must have length 1.")
        } else {
            return(v)
        }
    }
}

#' @rdname indepcopcond
#' @export
qcondindepcop <- function(tau,u,cpar=0){
    nu <- length(u)
    ntau <- length(tau)
    if (ntau == 1) {
        return(rep(tau, nu))
    } else {
        if (nu != ntau & nu != 1) {
            stop("u and tau must be the same length, or one must have length 1.")
        } else {
            return(tau)
        }
    }
}

#' @rdname indepcopcond
#' @export
qcondindep <- function(tau,u,cpar=0) qcondindepcop(tau,u,cpar=0)
