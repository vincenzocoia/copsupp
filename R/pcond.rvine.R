#' Conditional Distribution in a Regular Vine
#'
#' Evaluates the conditional distribution of a variable in a regular vine,
#' given values of the other variables.
#'
#' @param dat vector or matrix of observations (columns are variables).
#' @param cond Integer; the variable you wish to condition on (i.e. the
#' column number of \code{dat}).
#' @param A Vine array matrix, possibly truncated.
#' Variables should be labelled to correspond to
#' the order they appear in \code{dat}, not so that
#' they match the column number.
#' @param copmat Upper triangular matrix of copula names corresponding to
#' the edges in \code{A}.
#' @param cparmat Upper triangular matrix of copula parameters
#' corresponding to the copula families in \code{copmat}.
#' @param FXmarg List of (univariate) marginal cdfs of X_1, ..., X_p;
#' each should be vectorized. Or a single function if the cdf is all the same.
#' @param FYmarg Marginal cdf of Y, vectorized.
#' @details This function could do one of two things, depending on the
#' scenario.
#'
#' \itemize{
#'  \item If variable \code{cond} is not a leaf (that is, a vine array cannot
#'  be written with it at the end), then the conditional density is integrated.
#'  \item If variable \code{cond} is a leaf, then Bo's
#'  \code{rVineTruncCondCDF} function in the \code{copreg} package is
#'  used to compute the conditional cdf.
#' }
pcond.rvine <- function(dat, cond, A, copmat, cparmat, Fmarg = identity) {
    if (is.vector(dat)) dat <- matrix(dat, nrow = 1)
    p <- ncol(A)
    ntrunc <- nrow(A) - 1
    ## Uniformize data:
    if (length(Fmarg) == 1) {
        udat <- Fmarg(dat)
    } else {
        udat <- dat
        for (col in 1:ncol(dat)) udat[, col] <- Fmarg[[col]](dat[, col])
    }
    ## Is cond a leaf? If so, get the vine array with it as a leaf.
    Aleaf <- releaf.varray(A, leaf = cond)
    if (is.null(Aleaf)) {
        res <- apply(dat, 1, function(row) {
            dens_ <- function(xcond){
                x <- row
                x[cond] <- xcond
                dR(x, A, copmat, cparmat)
            }
            dens <- Vectorize(dens_)
            integrate(dens, 0, row[cond])$value / integrate(dens, 0, 1)$value
        })
    } else {
        ## We'll need to re-arrange the copmat and cparmat to match Aleaf.
        copmat <- reform.copmat(copmat, Aleaf, A)
        cparmat <- reform.copmat(cparmat, Aleaf, A)
        if (is.dvine(Aleaf)) {

        }
        ## Re-order data so that it's in order of the vine array, and re-label
        ##  the variables in the vine.
        perm <- diag(Aleaf)
        Aleafperm <- varrayperm(Aleaf, perm)
        datperm <- dat[, perm]
        ## Use Bo's function
        rVineTruncCondCDF
    }
    res

}
