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
#' @param A If you know the vine array that you want to use, put it here.
#' Labels should correspond to column numbers of \code{xdat}. Otherwise,
#' leave it \code{NULL}.
#' @param families A vector of copula family names to try
#' fitting (will also consider their rotations/reflections). Limited to
#' those families available in \code{VineCopula} package, listed in
#' \code{\link{BiCopSelect}}.
#' @param ... Other arguments to pass to \code{VineCopula::RVineCopSelect}.
#' @note Because of some conflict between the deprecated package \code{igraph0}
#' (used by package \code{CopulaModel})
#' and the newer package \code{igraph} (used by package \code{VineCopula}),
#' this function can only be run a maximum of one time.
#' Run it more and you'll get an error.
#' @return A list with three entries:
#'
#' \itemize{
#'      \item \code{$cdf}: List of length \code{ncol(xdat)} of the marginal
#'      distribution functions, as input in \code{margs}.
#'      \item \code{$A}: Vine array, truncated to \code{ntrunc}.
#'      \item \code{$copmat}: \code{ntrunc x ncol(A)} upper-triangular
#'      matrix of copula model names.
#'      \item \code{$cparmat}: \code{ntrunc x ncol(A)} upper-triangular
#'      matrix of copula parameters. If at least one copula family has more or
#'      less parameters than 1, each entry is a list of length one containing
#'      the vector of copula parameters for that copula family.
#' }
#' 
#' Note that the \code{$cdf} output is listed so that the position of the cdf
#' in the list corresponds to the integer labels in \code{$A}. 
#' This might change in the future if it's more convenient or sensible
#' to only listing the cdfs for variables in \code{$A}.
#' @details For the \code{familyset} argument, the default is almost all of
#' the families available. It just doesn't include the Tawn copula families.
#'
#' Note that you'll need the \code{igraph} package installed.
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
fit.rvine <- function(xdat, vars = 1:ncol(xdat), ntrunc = ncol(xdat)-1, 
                      margs = identity, A = NULL,
                      families = c("bvncop","bvtcop","mtcj","gum","frk","joe","bb1","bb7","bb8"), ...) {
    familyset = sort(unique(c(copname2num(families), recursive = TRUE)))
    if (is.vector(xdat)) xdat <- matrix(xdat, ncol = 1)
    if (is.data.frame(xdat)) xdat <- as.matrix(xdat)
    p_all <- ncol(xdat)
    n <- nrow(xdat)
    p <- length(vars)
    ## Marginals:
    if (length(margs) == 1) margs <- rep(list(margs), p_all)
    if (p == 1){
        return(list(cdf=margs, A=matrix(vars), copmat=matrix(""), cparmat=matrix(0)))
    }
    ## Uniformize and subset data
    for (col in vars) xdat[, col] <- margs[[col]](xdat[, col])
    xdat <- xdat[, vars]
    ## Get vine array to input to RVineCopSelect. Call it `A1` instead of `A`
    ##  because I'll use the output matrix of RVineCopSelect (which *should*
    ##  always be the same anyway).
    if (is.null(A)) {
        if (p == 2) {
            A1 <- matrix(c(1,0,1,2), ncol = 2) # Needs to have labels = 1:2.
        } else {
            ## Get correlation matrix
            cormat <- cor(qnorm(xdat))
            ## Choose vine array
            library(igraph)
            arrayfit <- CopulaModel::gausstrvine.mst(cormat, ntrunc)
            A1 <- arrayfit$RVM$VineA
        }
    } else {
        A1 <- A
    }
    
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
    if (!is.matrix(parmat1)) parmat1 <- matrix(parmat1, ncol = p)
    if (!is.matrix(parmat2)) parmat2 <- matrix(parmat2, ncol = p)
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
    ## It seems that, when RVineCopSelect fits a 90- or 270-degree rotated
    ##  copula, it also makes the parameters negative. Need to fix this.
    ## joe; mtcj; gum; bb1; bb6; bb7; bb8
    for (i in 1:nrow(copmat)) for (j in (i+1):ncol(copmat)) {
        thiscop <- copmat[i, j]
        ## Is this copula model 90- or 270-degree rotated?
        lastchar <- substring(thiscop, first=nchar(thiscop), last=nchar(thiscop))
        if (lastchar == "u" | lastchar == "v") {
            ## Copula name ends in "u" or "v". Could that mean that it's "rotated"?
            wouldbe_pcop <- paste0("p", substring(thiscop, first=1, last=nchar(thiscop)-1))
            if (exists(wouldbe_pcop)) rotated <- TRUE else rotated <- FALSE
            ## Now that we know if it's rotated, flip the parameters if need be.
            if (rotated) {
                if (is.list(cparmat[1,1])) {
                    cparmat[i, j][[1]] <- -cparmat[i, j][[1]]
                } else {
                    cparmat[i, j] <- -cparmat[i, j]
                }
            }
        }
    }
    ## Output results
    Avars <- varray.vars(A)
    list(cdf=margs,
         A=relabel.varray(A, vars[Avars]), 
         copmat=copmat, 
         cparmat=cparmat)
}
