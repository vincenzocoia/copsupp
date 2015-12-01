#' Fit a Vine Model
#'
#' Fits a joint distribution for the data using a vine copula. The vine array
#' is first chosen using the minimum spanning tree algorithm using
#' the function \code{gausstrvine.mst} in the \code{CopulaModel}
#' package, then the pairwise
#' copula models are chosen and fit using \code{RVineCopSelect} in the
#' \code{VineCopula} package.
#'
#' @param xdat Matrix of data; columns represent variables, and rows observations.
#' @param vars Vector of integers specifying the column numbers of \code{xdat}
#' to fit a model to. Default is all variables.
#' @param ntrunc Integer, either \code{1, 2, ...,ncol(xdat)-1},
#' of the truncation level of the vine to be fit.
#' @param margs List of vectorized functions of the univariate marginal cdf's of
#' the data, in the order of the columns. Or a single function if the cdf's
#' are all the same.
#' @param familyset A vector of integer codes of the copula families to try
#' fitting. See \code{VineCopula::RVineCopSelect} for a full list.
#' @param ... Other arguments to pass to \code{VineCopula::RVineCopSelect}.
#' @note Because of some conflict between the deprecated package \code{igraph0}
#' (used by package \code{CopulaModel})
#' and the newer package \code{igraph} (used by package \code{VineCopula}),
#' this function can only be run a maximum of one time.
#' Run it more and you'll get an error.
#' @return A list with three entries:
#'
#' \itemize{
#'      \item \code{$A}: Vine array, truncated to \code{ntrunc}.
#'      \item \code{$copmat}: \code{ntrunc x ncol(A)} upper-triangular
#'      matrix of copula model names.
#'      \item \code{$cparmat}: \code{ntrunc x ncol(A)} upper-triangular
#'      matrix of copula parameters. If at least one copula family has more or
#'      less parameters than 1, each entry is a list of length one containing
#'      the vector of copula parameters for that copula family.
#' }
#' @details For the \code{familyset} argument, the default is almost all of
#' the families available. It just doesn't include the Tawn copula families.
#'
#' Note that you'll need the \code{igraph0} package installed.
#' @examples
#' ## Get some simulated data:
#' set.seed(152)
#' ntrunc <- 2
#' p <- 4
#' A0 <- truncvarray(CopulaModel::Dvinearray(p), ntrunc)
#' copmat0 <- makeuppertri("frk", ntrunc, p, "")
#' cparmat0 <- makeuppertri(3, ntrunc, p)
#' dat <- fvinesim(100, A0, copmat0, cparmat0)
#'
#' ## Fit a model to the data:
#' fit.rvine(dat, ntrunc=ntrunc)
#' fit.rvine(dat, c(4, 2, 3))
#' @export
fit.rvine <- function(xdat, vars = 1:ncol(xdat), ntrunc = ncol(xdat)-1, margs = identity,
                 familyset = c(1:10,13,14,16:20,23,24,26:30,33,34,36:40), ...) {
    if (is.vector(xdat) | ncol(xdat) == 1){
        list(A=matrix(1), copmat=matrix(""), cparmat=matrix(0))
    }
    p_all <- ncol(xdat)
    n <- nrow(xdat)
    if (length(margs) == 1) margs <- rep(list(margs), p_all)
    ## Uniformize and subset data
    for (col in vars) xdat[, col] <- margs[[col]](xdat[, col])
    xdat <- xdat[, vars]
    p <- length(vars)
    ## Get correlation matrix
    cormat <- cor(qnorm(xdat))
    ## Choose vine array
    # library(igraph0)
    library(igraph)
    arrayfit <- CopulaModel::gausstrvine.mst(cormat, ntrunc)
    A1 <- arrayfit$RVM$VineA
    ## Now get and fit copulas
    capture.output(vinefit <- VineCopula::RVineCopSelect(xdat,
                                                         familyset = familyset,
                                                         Matrix = A1,
                                                         trunclevel = ntrunc,
                                                         ...))
    ## Extract things.
    #### Vine array -- it should be the same as A1. Truncate it if need be.
    A <- vinefit$Matrix[p:1, p:1]
    if (!identical(A, A1))
        warning(paste("The vine array output by RVineCopSelect is different",
                      "than what was input. Using the output anyway, but you",
                      "might want to investigate why."))
    A <- truncvarray(A, ntrunc)
    #### copmat
    copmatind <- vinefit$family[(p:1)[1:ntrunc], p:1]
    if (!is.matrix(copmatind)) copmatind <- matrix(copmatind, ncol = p)
    copmat <- apply(copmatind, 1:2, copnum2name)
    copmat[!upper.tri(copmat)] <- ""
    #### cparmat
    parmat1 <- vinefit$par[(p:1)[1:ntrunc], p:1]
    parmat2 <- vinefit$par2[(p:1)[1:ntrunc], p:1]
    parvec <- numeric(0)
    len <- integer(0)
    for (i in 1:ntrunc) for (j in (i+1):p) {
        if (parmat1[i, j] == 0) {
            len <- c(len, 0)
        } else {
            parvec <- c(parvec, parmat1[i, j])
            if (parmat2[i, j] != 0) {
                parvec <- c(parvec, parmat2[i, j])
                len <- c(len, 2)
            } else {
                len <- c(len, 1)
            }
        }
    }
    if (all(len == 1)) {
        cparmat <- makeuppertri(parvec, nrow = ntrunc, ncol = p)
    } else {
        cparmat <- makeuppertri.list(parvec, len, nrow = ntrunc, ncol = p)
    }
    ## Output results
    Avars <- varray.vars(A)
    list(A=relabel.varray(A, vars[Avars]), copmat=copmat, cparmat=cparmat)
}
