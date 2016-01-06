#' Fit copula models to a vine array.
#'
#' Choose and fit copula models on a pre-specified vine array. Intended for
#' internal use, as opposed to \code{RVineCopSelect}. Currently still
#' uses \code{VineCopula} through \code{BiCopSelect}.
#'
#' @param dat Data matrix with Uniform margins.
#' @param A Vine array matrix, possibly truncated.
#' @param copmat Pre-specified copula families in the form of an upper-triangular
#' matrix. Put \code{NA} to leave the edge unspecified. \code{NULL} for
#' fully unspecified.
#' @param cparmat Pre-specified copula parameters corresponding to some of the
#' specified copulas in \code{copmat}. Put \code{NA} in place of parameters to
#' leave them unspecified. \code{NULL} for fully unspecified.
#' @param famililes Vector of candidate copula family names.
#' @import VineCopula
#' @note Expecting smart input. So, don't input a vine array that has no edges,
#' and don't specify a parameter for an unspecified copula family.
fitvarray <- function(dat, A, copmat=NULL, cparmat=NULL,
                      families = c("indepcop", "bvncop","bvtcop","mtcj","gum",
                                   "frk","joe","bb1","bb7","bb8")) {
    ntrunc <- nrow(A) - 1
    d <- ncol(A)
    if (is.null(copmat))
        copmat <- makeuppertri(NA, ntrunc, d)
    if (is.null(cparmat))
        cparmat <- makeuppertri.list(NA, rep(1, ntrunc*d - choose(ntrunc+1,2)), ntrunc, d)
    startparmat <- cparmat
    ## Number of parameters for these copula families:
    allfams <- unique(c(copmat[upper.tri(copmat)], families))
    numpars <- sapply(allfams, function(cop_) length(cparspace(cop_, FALSE)$lower))
    names(numpars) <- allfams
    ## Go through edges and choose the best copula models. Do so by adding
    ##  variables one at a time.
    for (j in 2:d) for (i in 1:(j-1)) {
        ## Get conditional set:
        condset <- A[seq_len(i-1), j]
        ## Pairs:
        V1 <- A[i, j]
        V2 <- A[j, j]
        ## Get copula families
        if (is.na(copmat[i, j])) thesefams <- families else thesefams <- copmat[i, j]
        ## Fit model(s) if parameter(s) not specified.
        if (any(is.na(cparmat[i, j][[1]]))) for (cop_ in thesefams) {
            ## Get logdcop function
            thislogdcop <- get(paste0("logd", cop_))
            ## Get those parts of the parameter that are specified:
            thisnumpar <- numpars[cop_]

            obj <- function(cpar) -sum(logdcop(dat[, V1], dat[, V2]))
        }



    }
}
