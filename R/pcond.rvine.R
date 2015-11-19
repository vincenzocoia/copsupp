#' Conditional Distribution in a Regular Vine
#'
#' Evaluates the conditional distribution of a variable in a regular vine,
#' given values of the other variables. WARNING: Algorithm for non-d-vine
#' may not be working properly.
#'
#' @param dat vector or matrix of observations (columns are variables).
#' @param cond Integer; the variable you wish to condition on (i.e. the
#' column number of \code{dat}).
#' @param A Vine array matrix, possibly truncated.
#' Variables should be labelled to correspond to
#' the order they appear in \code{dat}, not so that
#' they match the column number of \code{A}.
#' @param copmat Upper triangular matrix of copula names corresponding to
#' the edges in \code{A}.
#' @param cparmat Upper triangular matrix of copula parameters
#' corresponding to the copula families in \code{copmat}.
#' @param Fmarg List of (univariate) marginal cdfs corresponding to
#' the columns in \code{dat};
#' each should be vectorized. Or a single function if the cdf is all the same.
#' @param .print Logical; should the function output how it goes about
#' finding the conditional distribution?
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
#' @return A vector of length = the number of observations in \code{dat},
#' representing the evaluated conditional distribution of variable \code{cond}
#' given the other variables in \code{A}.
#' @examples
#' ## D-Vine example
#' A <- CopulaModel::Dvinearray(5)
#' A <- relabel.varray(A, c(1, 5, 4, 3, 2))
#' A <- truncvarray(A, 2)
#' copmat <- makeuppertri("bvncop", 2, 5)
#' cparmat <- makeuppertri(c(1:7/10), 2, 5, byRow = FALSE)
#' udat <- fvinesim(10, A, copmat, cparmat)
#' pcond.rvine(udat, 5, A, copmat, cparmat, .print=T)  # integrates vine density.
#' pcond.rvine(udat, 2, A, copmat, cparmat, .print=T)  # computes from D-vine formula
#'
#' ## C-Vine example
#' A <- CopulaModel::Cvinearray(5)
#' A <- truncvarray(A, 2)
#' udat <- fvinesim(10, A, copmat, cparmat)
#' pcond.rvine(udat, 3, A, copmat, cparmat, .print=T)  # computes from general R-vine algo
#'
#' ## Array doesn't have to involve all data:
#' A <- CopulaModel::Dvinearray(5)
#' A <- truncvarray(A, 2)
#' A <- rvinesubset(A, 3:5)
#' copmat <- makeuppertri("frk", 2, 3, "")
#' cparmat <- makeuppertri(3:1, 2, 3)
#' pcond.rvine(1:5/10, 3, A, copmat, cparmat)
#' pcond.rvine(1:5/10, 4, A, copmat, cparmat)
#' ## are the same as...
#' A <- CopulaModel::Dvinearray(3)
#' A <- truncvarray(A, 2)
#' pcond.rvine(3:5/10, 1, A, copmat, cparmat)
#' pcond.rvine(3:5/10, 2, A, copmat, cparmat)
#' @export
pcond.rvine <- function(dat, cond, A, copmat, cparmat, Fmarg = identity,
                        .print = FALSE) {
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
        if (.print) print(paste0("cond=", cond, " is not a leaf. ",
                                 "Obtaining conditional distribution by integration."))
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
        ## We'll need to re-arrange the data so that it's in order of the
        ##  vine array.
        vars <- varray.vars(Aleaf)
        udat <- udat[, vars]
        if (!is.matrix(udat)) udat <- matrix(udat, ncol = length(vars))
        if (ncol(udat) <= 2) {
            ## NOTE: The pcondD (and qcondD) functionality requires the copmat
            ##  and cparmat to be in natural order. I won't bother with it,
            ##  but since rVineTruncCondCDF() won't accept a vine with 2
            ##  variables, and pcondD will, I'll use pcondD for that case.
            if (.print) print(paste0("using `pcondD.generic()`. because there",
                                     " are two variables left."))
            res <- apply(udat, 1, function(row){
                pcondD.generic(row, p, -(p-1),
                               copmat = copmat, cparmat = cparmat)
            })
        } else {
            if (.print) print(paste0("cond=", cond, " is a leaf, not of a D-vine. ",
                                     "using `copreg::rVineTruncCondCDF()`."))
            ## Use Bo's function
            ## Re-label the variables in the vine so that they're 1:p, in that
            ##  order.
            Aleaf <- relabel.varray(Aleaf)
            ## Fill-in vine array so it's p x p
            Aleaf <- rbind(Aleaf, matrix(0, nrow = p - ntrunc - 1, ncol = p))
            diag(Aleaf) <- 1:p
            ## parvec
            # parvec <- c(t(cparmat), recursive = TRUE)
            parvec <- t(cparmat)[lower.tri(t(cparmat))]
            ## pcondmat
            pcondmat <- apply(copmat, 1:2, function(cop) paste0("pcond", cop))
            pcondmat[!upper.tri(pcondmat)] <- ""
            ## np
            if (is.list(cparmat[1,1])) {
                np <- apply(cparmat, 1:2, function(cpar) length(cpar[[1]]))
            } else {
                np <- makeuppertri(1, nrow = nrow(cparmat), ncol = ncol(cparmat))
            }
            ## Call:
            library(CopulaModel)
            res <- copreg::rVineTruncCondCDF(parvec = parvec,
                                             udat = udat,
                                             A = Aleaf,
                                             ntrunc = ntrunc,
                                             pcondmat = pcondmat,
                                             np = np)
       }
    }
    res
}
