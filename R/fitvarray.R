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

#' @param edges Vector of c(node, connection1, connection2, etc.)
fitlayer <- function(dat, basevine, edges, cops = NULL, cpars = NULL,
                     families = c("indepcop", "bvncop","bvtcop","mtcj","gum",
                                  "frk","joe","bb1","bb7","bb8")) {
    d <- length(edges)
    if (is.null(cops)) cops <- rep(NA, d)
    if (is.null(cpars)) cpars <- rep(list(NA), d)
    if (!is.list(cops)) cops <- as.list(cops)
    if (!is.list(cpars)) cpars <- as.list(cpars)
    startpars <- cpars
    ## Number of parameters for these copula families:
    allfams <- unique(c(copmat[upper.tri(copmat)], families))
    numpars <- sapply(allfams, function(cop_) length(cparspace(cop_, FALSE)$lower))
    names(numpars) <- allfams
    ## Go through edges and choose the best copula models. Do so by adding
    ##  variables one at a time.
    # for (j in 2:d) for (i in 1:(j-1)) {
    u <- dat[, edges[2]]
    v <- dat[, edges[1]]
    for (i in 2:length(edges) - 1) {

#         ## Pairs:
#         V1 <- A[i, j]
#         V2 <- A[j, j]
        ## Get candidate copula families
        if (any(is.na(cops[[i]]))) cand <- families else cand <- cops[[i]]
        ## Fit model(s) if parameter(s) not specified.
        thiscpars <- cpars[[i]]
        if (any(is.na(thiscpars))) {
            ## Choose copula and get starting values
            startfit <- VineCopula::BiCopSelect(u, v, familyset = copname2num(cops[[i]]))
            vinecopcpars <- c(startfit$par, startfit$par2)
            if (vinecopcpars[2] == 0) vinecopcpars <- vinecopcpars[-2]
            bestcop <- copnum2name(startfit$family)
            ## Which parameters are pre-specified?
            pre_spec <- !is.na(thiscpars)
            if (any(pre_spec)) {  # If any part of the parameter is pre-specified, re-do fitting.
                ## Get objective function
                logdcop <- get(paste0("logd", bestcop))
                obj <- function(theta) {
                    fulltheta <- vinecopcpars # Dummy
                    fulltheta[!pre_spec] <- theta
                    fulltheta[pre_spec] <- thiscpars[pre_spec]
                    -sum(logdcop(u, v, fulltheta))
                }
                finalfit <- nlm(obj, vinecopcpars[!pre_spec])
                thiscpars[!pre_spec] <- finalfit$estimate
            } else {
                thiscpars <- vinecopcpars
            }
        }
        ## This edge is now fit. Update and get new (u,v) variables if there's still more.



    }
}
