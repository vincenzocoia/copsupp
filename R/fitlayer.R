#' Fit a new layer
#'
#' Choose and fit copula models on a new layer. The edge ("array column") must be
#' pre-specified. Intended for internal use.
#'
#' @param dat Data matrix with Uniform margins.
#' @param basevine Object of type "rvine" of the already-fit base vine for which
#' the new layer is to be applied.
#' @param edges Vector -- new column of vine array (with node appearing first)
#' @param cops Vector or list of pre-specified
#' copula families for each edge. Put \code{NA}
#' to leave the edge unspecified. \code{NULL} for
#' fully unspecified. You're allowed to put more than one family
#' as candidates.
#' @param cpars Pre-specified copula parameters corresponding to some of the
#' specified copulas in \code{cops}. Put \code{NA} in place of parameters to
#' leave them unspecified. \code{NULL} for fully unspecified.
#' @param families Vector of candidate copula family names for those that are
#' \code{NA} or \code{NULL}.
#' @import VineCopula
#' @note Expecting smart input. So, ensure that \code{edges} has length at
#' least 2, and that edges[-1] are variables in \code{basevine}, and that
#' \code{cpars} are only specified when there's only one copula family to
#' choose from.
#' @details Edges are fit so that edges[1] is the "V" variable. So, copulas
#' are fit to (edges[2], edges[1]), then (edges[3],
#' edges[1]) | edges[2], etc. That's because when computing edges[1]|others,
#' "pcond" can be used instead of "pcond12".
#' @return List of fitted \code{$cops} and \code{$cpars}.
#' @export
fitlayer <- function(dat, basevine, edges, cops = NULL, cpars = NULL,
                     families = c("indepcop", "bvncop","bvtcop","mtcj","gum",
                                  "frk","joe","bb1","bb7","bb8")) {
    d <- length(edges)
    if (is.null(cops)) cops <- rep(NA, d-1)
    if (is.null(cpars)) cpars <- rep(list(NA), d-1)
    if (!is.list(cops)) cops <- as.list(cops)
    if (!is.list(cpars)) cpars <- as.list(cpars)
    startpars <- cpars
    ## Go through edges and choose the best copula models. Do so by adding
    ##  variables one at a time.
    # for (j in 2:d) for (i in 1:(j-1)) {
    u <- dat[, edges[2]]
    v <- dat[, edges[1]]
    for (i in 2:length(edges) - 1) {
        ## Get candidate copula families
        if (any(is.na(cops[[i]]))) cand <- families else cand <- cops[[i]]
        ## Fit edge if cpar is not fully specified.
        if (all(is.na(cpars[[i]]))) {
            thisfit <- fitbicop_lh(u, v, families = cand)
        } else {
            thisfit <- fitbicop_lh(u, v, families = cand, cpar = cpars[[i]])
        }
        ## This edge is now fit. Update and get new (u,v) variables if there's still more.
        cops[[i]] <- thisfit$cop
        cpars[[i]] <- thisfit$cpar
        if (i+1 < length(edges)) {
            pcond <- get(paste0("pcond", thisfit$cop))
            v <- pcond(v, u, cpars[[i]])
            u <- pcondrvine(dat, basevine, cond = edges[i+2], vbls = edges[2:(i+2)])
        }
    }
    ## Output fits
    list(cops = c(cops, recursive = TRUE), cpars = cpars)
}

#' Temporary Title
#'
#' @param cops List with vector entries of the copula families to try fitting
#' for each edge.
#' @param cpars List with entries being the list of parameters you want to enforce
#' on the copula models. \code{NULL} for no restrictions; \code{NULL} entry
#' for no restrictions on that edge; \code{NA} in a parameter vector
#' means that parameter needs fitting.
#' @param QY Vectorized function of the fitted quantile function of Y.
#' @param families Vector of copula families.
#' For those edges in \code{cop} that don't have copula families
#' specially selected (i.e. have \code{NULL} entries), these families are used.
fitlayer_cnqr <- function(dats, edges, basevine, tauset, QY, stoponbest = TRUE,
                          w = function(x) 1, cops = NULL, cpars = NULL,
                          families = c("indepcop", "bvncop","bvtcop","mtcj","gum",
                                       "frk","joe","bb1","bb7","bb8")) {
    ## 1. Modify and get info from input
    dattr <- dats$tr
    datval <- dats$val
    ## Fill-in cops if not done already.
    nfam <- length(families)
    nedge <- length(edges) - 1  # Number of edges
    if (is.null(cpar)) cpar <- rep(list(NULL), nfam)
    cpar <- lapply(cpar, function(fams){
        if (is.null(fams)) families else fams
    })
    ## Make cpars a list if it's not already.
    if (!is.list(cpars)) {
        cpars <- rep(list(NULL), nedge)
    }
    ## 2. Get sequential predictor cdfs
    ucondtr <- matrix(dattr[, edges[2]], ncol = 1)
    ucondval <- matrix(datval[, edges[2]], ncol = 1)
    for (i in 2+seq_len(nedge - 1)) {
        ## Training
        nextu <- pcondrvine(dattr, basevine,
                            cond = edges[i],
                            vbls = edges[2:i])
        ucondtr <- cbind(ucondtr, nextu)
        ## Validation
        nextu <- pcondrvine(datval, basevine,
                            cond = edges[i],
                            vbls = edges[2:i])
        ucondval <- cbind(ucondval, nextu)
    }
    ## 3. Fit each edge, one at a time.
    for (i in 1:nedge) {
        ## BEGIN i'th EDGE
    }




    ## Initiate output
    valscore <- numeric(0)
    fittedcpars <- list()
    fittedcops <- character()  # We need this because the copula could be flipped.
    ## Fit each model
    for (i in 1:nfam) {
        cop <- families[i]
        ## Do we need to estimate any parameters? No, if the parameters are already
        ##  fully specified, or if there are no parameters.
        pre_cpar <- cpar[[i]]
        do_estimate <- TRUE
        if (!is.null(pre_cpar)) if (all(!is.na(pre_cpar))) {
            do_estimate <- FALSE
        }
        if (cop == "indepcop") {
            do_estimate <- FALSE
            pre_cpar <- numeric(0)
        }
        ## Were the parameters only partially specified?
        ## (CHECK)
        ## Estimate if need be.
        if (!do_estimate) {
            ## No need to estimate. Just get forecasts on the validation set.
            qcond <- get(paste0("qcond", cop))
            vhat_val <- sapply(taus, function(tau_) qcond(tau_, uval, pre_cpar))
        } else {
            ## We'll have to estimate. VineCopila::BiCopSelect has smart starting values,
            ##  so I'll use those to get starting values.
            startfit <- BiCopSelect(u, v, familyset = copname2num(cop)[[1]])
            ## Get proper copula rotation
            cop <- copnum2name(startfit$family)
            qcond <- get(paste0("qcond", cop))
            ## Extract starting values
            cparstart <- c(startfit$par, startfit$par2)
            if (cparstart[2] == 0) cparstart <- cparstart[1]
            ## Get forecasting family:
            ## ***NEED TO PUT IN PRE-SPECIFIED PORTION OF PARAMETER.
            forecast_fam <- function(u, theta) {
                sapply(taus, function(tau_) qcond(tau_, u, theta))
            }
            ## Choose theta to minimize score on training set:
            vhat_train <- function(theta) forecast_fam(u, theta)
            obj <- function(theta) scoreq(v, vhat_train(theta), taus)
            cparhat <- rnlm(obj, cparstart, cparspace(cop))$estimate
            ##
            ## Get forecast on validation set:
            vhat_val <- forecast_fam(uval, cparhat)
        }
        ## Get score from validation set:
        thisscore <- scoreq(vval, vhat_val, taus)
        ## Update outputs.
        valscore[i] <- thisscore
        fittedcpars[[i]] <- cparhat
        fittedcops[i] <- cop
    }
    ## Output
    list(cops = fittedcops, cpars = fittedcpars, scorehat = valscore)
}
