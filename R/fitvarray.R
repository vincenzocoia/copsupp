# #' Fit copula models to a vine array.
# #'
# #' Choose and fit copula models on a pre-specified vine array. Intended for
# #' internal use, as opposed to \code{RVineCopSelect}. Currently still
# #' uses \code{VineCopula} through \code{BiCopSelect}.
# #'
# #' @param dat Data matrix with Uniform margins.
# #' @param A Vine array matrix, possibly truncated.
# #' @param copmat Pre-specified copula families in the form of an upper-triangular
# #' matrix. Put \code{NA} to leave the edge unspecified. \code{NULL} for
# #' fully unspecified.
# #' @param cparmat Pre-specified copula parameters corresponding to some of the
# #' specified copulas in \code{copmat}. Put \code{NA} in place of parameters to
# #' leave them unspecified. \code{NULL} for fully unspecified.
# #' @param famililes Vector of candidate copula family names.
# #' @import VineCopula
# #' @note Expecting smart input. So, don't input a vine array that has no edges,
# #' and don't specify a parameter for an unspecified copula family.
# fitvarray <- function(dat, A, copmat=NULL, cparmat=NULL,
#                       families = c("indepcop", "bvncop","bvtcop","mtcj","gum",
#                                    "frk","joe","bb1","bb7","bb8")) {
#     ntrunc <- nrow(A) - 1
#     d <- ncol(A)
#     if (is.null(copmat))
#         copmat <- makeuppertri(NA, ntrunc, d)
#     if (is.null(cparmat))
#         cparmat <- makeuppertri.list(NA, rep(1, ntrunc*d - choose(ntrunc+1,2)), ntrunc, d)
#     startparmat <- cparmat
#     ## Number of parameters for these copula families:
#     allfams <- unique(c(copmat[upper.tri(copmat)], families))
#     numpars <- sapply(allfams, function(cop_) length(cparspace(cop_, FALSE)$lower))
#     names(numpars) <- allfams
#     ## Go through edges and choose the best copula models. Do so by adding
#     ##  variables one at a time.
#     for (j in 2:d) for (i in 1:(j-1)) {
#         ## Get conditional set:
#         condset <- A[seq_len(i-1), j]
#         ## Pairs:
#         V1 <- A[i, j]
#         V2 <- A[j, j]
#         ## Get copula families
#         if (is.na(copmat[i, j])) thesefams <- families else thesefams <- copmat[i, j]
#         ## Fit model(s) if parameter(s) not specified.
#         if (any(is.na(cparmat[i, j][[1]]))) for (cop_ in thesefams) {
#             ## Get logdcop function
#             thislogdcop <- get(paste0("logd", cop_))
#             ## Get those parts of the parameter that are specified:
#             thisnumpar <- numpars[cop_]
#
#             obj <- function(cpar) -sum(logdcop(dat[, V1], dat[, V2]))
#         }
#
#
#
#     }
# }

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
