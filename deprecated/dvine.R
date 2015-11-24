#' Conditional D-Vine Distribution: Generic variables
#'
#' For the vine X_1--X_2--X_3--..., these functions return the joint cdf or
#' quantile function of X_i given the \code{num} variables either above or
#' below X_i in the vine. Meant primarily for internal use.
#'
#' @param x vector containing at least the values of rv's \eqn{X_i} to \eqn{X_{i+num}},
#' where the position in the vector matches the variable index (i.e. i,...,i+n).
#' @param i integer index of unconditioned variable (\eqn{X_i}).
#' @param num integer; number of variables directly upstream or downstream
#' of X_i to condition on. For upstream (i+1, i+2, ...), make
#' \code{sign(num)>0}; for downstream (i-1, i-2, ...), make \code{sign(num)<0}.
#' @param copmat Upper triangular \code{(num trees)x(i+num)} matrix of
#' copula names (like "frk" and "gum") as strings. Must be listed in "reading
#' order" (for example, natural order won't work).
#' @param cparmat Upper triangular matrix of copula parameters corresponding
#' to the copula families in \code{copmat}.
#' @param Fmarg List of marginal distributions (univariate functions),
#' whose position in the list corresponds to the variable index.
#' Must contain at least those marginals i:(i+n). Could be a single such function
#' if it's the same for all data.
#' @return A single numeric value representing the evaluated conditional
#' cdf or quantile function.
#' @note You might want to load the CopulaModel package to load the
#' copula families.
#' @details The copulas link (downstream variable, upstream variable), not
#' the other way around. This doesn't make a difference if the copula is "U-V"
#' symmetric.
#' @examples
#' ## Joint distribution info
#' library(CopulaModel)
#' copmat <- makeuppertri("frk", 5, 5)
#' cparmat <- makeuppertri(1:10, 5, 5)
#' Fmarg <- list(pnorm, pnorm, pnorm, pnorm, pnorm)
#'
#' ## cdf of X1|X2,...,X5 at -2|-1:2, with 2-truncation:
#' pcondD.generic(-2:2, 1, 4, copmat[1:2, ], cparmat[1:2, ], Fmarg)
#'
#' ## cdf of X4|X2,X3 at 0|1,1:
#' pcondD.generic(c(NA, 1, 1, 0), 4, -2, copmat, cparmat, Fmarg) # X_1 value doesn't matter.
#' @rdname dvine.generic
#' @export
pcondD.generic <- function(x, i, num, copmat, cparmat, Fmarg = identity){
    if (length(Fmarg) == 1) Fmarg <- rep(list(Fmarg), max(i, i+num))
    if (num == 0) {
        ## There's no conditioned variable, so just call on the marginal.
        Fmarg[[i]](x[i])
    } else { ## There are conditioned variables.
        ## One less variable to condition on in the next iteration:
        newnum <- sign(num)*(abs(num) - 1)
        ## Grab the copula
        minind <- min(i, i + num)
        maxind <- max(i, i + num)
        tree <- maxind - minind
        if (tree > nrow(copmat)) {
            thispcond <- pcondindepcop
            thiscpar <- integer(0)
        } else {
            thispcond <- paste0("pcond", copmat[tree, maxind])
            if (num > 0) {
                thispcond12 <- paste0(thispcond, "12")
                if (exists(thispcond12)) thispcond <- thispcond12
            }
            thispcond <- get(thispcond)
            thiscpar <- cparmat[tree, maxind]
            if (is.list(thiscpar)) thiscpar <- thiscpar[[1]]
        }
        ## Recursion formula
        thispcond(pcondD.generic(x, i,        newnum, copmat, cparmat, Fmarg),
                  pcondD.generic(x, i + num, -newnum, copmat, cparmat, Fmarg),
                  thiscpar)
    }
}

#' Conditional D-Vine Distribution
#'
#' For the vine Y--X_1--...--X_p,
#' these functions return the joint cdf or
#' quantile function of Y|X_1,...,X_p.
#'
#' @param y response(s) to evaluate the cdf at. Could be a vector (or matrix
#' with one column) to evaluate at multiple cdf's corresponding to the rows
#' of \code{x}.
#' @param x vector of covariates X_1,...,X_p; or, a matrix with different
#' covariate values (p columns).
#' @param copmat Upper triangular \code{ntrunc x p} matrix of
#' copula names (like "frk" and "gum") as strings in the D-vine array.
#' Must be listed in "reading
#' order" (for example, natural order won't work).
#' @param cparmat Upper triangular matrix of copula parameters corresponding
#' to the copula families in \code{copmat}.
#' @param FXmarg List of (univariate) marginal cdfs of X_1, ..., X_p;
#' each should be vectorized. Or a single function if the cdf is all the same.
#' @param FYmarg Marginal cdf of Y, vectorsized.
#' @note You might want to load the CopulaModel package to load the
#' copula families.
#' @details The copulas link (downstream variable, upstream variable), not
#' the other way around. This doesn't make a difference if the copula is "U-V"
#' symmetric.
#' @return Vector of numerics representing the evaluated
#' conditional cdfs or quantile functions.
#' @examples
#' ## Joint distribution info
#' library(CopulaModel)
#' copmat <- makeuppertri("frk", 5, 5)
#' cparmat <- makeuppertri(1:10, 5, 5)
#'
#' ## cdf of X1|X2,...,X5 at -2|-1:2
#' pcondD(-2, -1:2, copmat, cparmat, pnorm)
#' @rdname pcondD
#' @export
pcondD <- function(y, x, copmat, cparmat, FXmarg = identity, FYmarg = identity) {
    if (is.vector(x)) x <- matrix(x, nrow = 1)
    num <- ncol(x)
    if (length(FXmarg) == 1) FXmarg <- rep(list(FXmarg), num)
    x <- cbind(y, x)
    Fmarg <- c(list(FYmarg), FXmarg)
    apply(x, 1, function(row){
        pcondD.generic(row, 1, num, copmat, cparmat, Fmarg)
    })
}

#' @param tau The quantile index
#' @param Qi_marg Marginal quantile function of X_i (univariate function)
#' @examples
#' ## 0.6-quantile of X1|(X2,...,X5)=(2,2,2,2).
#' qcondD.generic(0.6, c(NA, 2, 2, 2, 2), 1, 4, copmat, cparmat, Fmarg, qnorm)
#'
#' ## 0.5-quantile of X2|(X3,X4)=(1,1). Notice only the 3rd and 4th positions
#' ##  in the x vector matters.
#' qcondD.generic(0.5, c(NA, 2, 1, 1, 99), 2, 2, copmat, cparmat, Fmarg, qnorm)
#' @rdname dvine.generic
#' @export
qcondD.generic <- function(tau, x, i, num, copmat, cparmat, Fmarg, Qi_marg){
    if (num == 0) {
        ## No conditioned variable, so take qf of X_i:
        Qi_marg(tau)
    } else { ## There are conditioned variables.
        ## One less variable to condition on in the next iteration:
        newnum <- sign(num)*(abs(num) - 1)
        ## Grab the copula
        maxind <- max(i, i + num)
        tree <- abs(num)
        thisqcond <- paste0("qcond", copmat[tree, maxind])
        if (num > 0) {
            thisqcond12 <- paste0(thisqcond, "12")
            if (exists(thisqcond12)) thispcond <- thispcond12
        }
        thisqcond <- get(thisqcond)
        thiscpar <- cparmat[tree, maxind]
        if (is.list(thiscpar)) thiscpar <- thiscpar[[1]]
        ## Recursion formula
        qcondD.generic(thisqcond(tau,
                                 pcondD.generic(x, i + num, -newnum, copmat, cparmat, Fmarg),
                                 thiscpar),
                       x, i, newnum, copmat, cparmat, Fmarg, Qi_marg)
    }
}

#' @param tau Vector of quantile indices to evaluate conditional quantile
#' functions
#' @param QY Quantile function of Y, vectorized.
#' @rdname pcondD
#' @export
qcondD <- function(tau, x, copmat, cparmat, FXmarg = identity, QY = identity) {
    if (is.vector(x)) x <- matrix(x, nrow = 1)
    p <- ncol(x)
    if (length(FXmarg) == 1) FXmarg <- rep(list(FXmarg), p)
    x <- cbind(NA, x)
    Fmarg <- c(NA, FXmarg)
    t(apply(x, 1, function(row){
        qcondD.generic(tau, row, 1, p, copmat = copmat, cparmat = cparmat,
                       Fmarg = Fmarg, Qi_marg = QY)
    }))
}
