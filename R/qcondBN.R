#' Conditional Distribution: Bayesian Network Method
#'
#' Finds the quantile function of a response \code{Y} given predictors
#' \code{X_1,...,X_p} through a Bayesian Network. The bivariate distributions
#' \code{(X_i, Y)|X_1,...,X_{i-1}} are
#' modelled using bivariate copula
#' models. The distribution of the response is assumed known/fitted, along
#' with the joint distribution of the predictors.
#'
#' @param tau Vector of quantile indices to evaluate at.
#' @param x Vector of predictor values \code{X_1, ..., X_p}. Length should
#' equal the number of variables conditioning on.
#' @param p Integer; number of variables \code{X_1,...,X_p} to condition on.
#' @param cop Vector of copula names (like "frk"), whose \code{i}th entry
#' links \code{(X_i, Y)|X_1,...,X_{i-1}}.
#' @param cpar List of copula parameters corresponding to \code{cop}. Or,
#' vector of parameters if each copula family in \code{cop} has one parameter.
#' @param fX Function that accepts a vector argument; the
#' density of \code{(X_1,...,X_p)}.
#' @param Fcond If you have the conditional cdf
#' \code{X_i|X_1,...,X_{i-1}} for some \code{i}'s, insert them here
#' in a list (in their respective components). Put \code{NULL}
#' for entries that you don't have the form for. You could also already
#' evaluate this function at \code{x} if you'd like.
#' @param QY Quantile function for the response.
#' @details If \code{Fcond} is not missing any cdf's, then there's no
#' need to specify \code{fX}.
#' @return Vector of conditional quantile functions evaluated at the quantile
#' indices \code{tau}.
#' @rdname qcondBN
#' @export
qcondBN <- function(tau, x, p, cop, cpar, fX, Fcond = NULL, QY) {
    if (p == 0) return(QY(tau))
    ## Get conditional cdf of the predictors, evaluated.
    if (is.null(Fcond[[p]])) {
        fp <- function(x_) {
            xvar <- x[1:p]
            xvar[p] <- x_
            fX(xvar)
        }
        fp <- Vectorize(fp)
        integ <- integrate(fp, -Inf, x[p])
        area <- integrate(fp, -Inf, Inf)
        if (integ$message != "OK")
            stop (paste0("Integrate error, variable ", p, ": ", integ$message))
        thisFcond <- integ$value / area$value
    } else {
        if (is.function(Fcond[[p]])) {
            if (p == 1) {
                thisFcond <- Fcond[[p]](x[p])
            } else {
                thisFcond <- Fcond[[p]](x[p], x[-p])
            }
        } else {
            thisFcond <- Fcond[[p]]
        }

    }
    ## Remove end variable from fX.
    fXold <- fX
    fX <- function(xlessp) {
        g <- function(t) fXold(c(xlessp, t))
        g <- Vectorize(g)
        integrate(g, -Inf, Inf)$value
    }
    ## Recursive form of conditional qf of response:
    thiscop <- get(paste0("qcond", cop[p]))
    if (is.list(cpar)) {
        thiscpar <- cpar[[p]]
    } else {
        thiscpar <- cpar[p]
    }
    thiscopfit <- function(arg1, arg2) thiscop(arg1, arg2, thiscpar)
    newtau <- thiscopfit(tau, thisFcond)
    return(qcondBN(newtau, x, p-1, cop, cpar, fX, Fcond, QY))
}

#' @rdname qcondBN
#' @export
qcondBNu <- function(tau, u, p, cop, cpar, fU, Fcond = NULL, QY) {
    if (p == 0) return(QY(tau))
    ## Get conditional cdf of the predictors, evaluated.
    if (is.null(Fcond[[p]])) {
        fp <- function(x_) {
            xvar <- u[1:p]
            xvar[p] <- x_
            fU(xvar)
        }
        fp <- Vectorize(fp)
        integ <- integrate(fp, 0, u[p])
        area <- integrate(fp, 0, 1)
        if (integ$message != "OK")
            stop (paste0("Integrate error, variable ", p, ": ", integ$message))
        thisFcond <- integ$value / area$value
    } else {
        if (is.function(Fcond[[p]])) {
            if (p == 1) {
                thisFcond <- Fcond[[p]](u[p])
            } else {
                thisFcond <- Fcond[[p]](u[p], u[-p])
            }
        } else {
            thisFcond <- Fcond[[p]]
        }

    }
    ## Remove end variable from fU.
    if (any(sapply(Fcond, is.null))) {
        fXold <- fU
        fU <- function(xlessp) {
            g <- function(t) fXold(c(xlessp, t))
            g <- Vectorize(g)
            integrate(g, 0, 1)$value
        }
    }
    ## Recursive form of conditional qf of response:
    thiscop <- get(paste0("qcond", cop[p]))
    if (is.list(cpar)) {
        thiscpar <- cpar[[p]]
    } else {
        thiscpar <- cpar[p]
    }
    thiscopfit <- function(arg1, arg2) thiscop(arg1, arg2, thiscpar)
    newtau <- thiscopfit(tau, thisFcond)
    return(qcondBNu(newtau, u, p-1, cop, cpar, fU, Fcond, QY))
}
