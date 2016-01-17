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

