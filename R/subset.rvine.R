#' Subset a Regular Vine
#'
#' Find the vine connecting a subset of variables from a bigger vine, if
#' it exists.
#'
#' @param rv A regular vine object.
#' @param select Vector of variable indices in \code{diag(A)} to subset,
#' if possible. The order of the variables does not matter.
#' @details Just a technicality:
#' by saying a subset "doesn't have an existing vine", I mean that
#' a vine can't be formed using nodes and edges from
#' the original -- not that the
#' joint distribution of the selected variables can't be created from a vine
#' (so as to say, for example, that the simplifying assumption of vines
#' doesn't hold for this distribution).
#' @return Returns a regular vine of the subsetted variables, with variables
#' ordered according to their order in \code{A}; or \code{NULL} if
#' the subset does not form a vine.
#' @examples
#' ## Setup a vine.
#' A <- CopulaModel::Dvinearray(5)
#' copmat <- makeuppertri(c("gum", "mtcj", "gal", "joe",
#'                          "frk", "gum", "bb7",
#'                          "bb1", "indepcop",
#'                          "bb8"), 4, 5, "")
#' cparmat <- makeuppertri.list(c(3, 2.5, 2, 1.5,
#'                                1, 1.3, 2, 2,
#'                                3, 4,
#'                                5, 0.5),
#'                                len = c(1,1,1,1,1,1,2,2,0,2),
#'                                4, 5)
#' (rv <- rvine(A, copmat, cparmat,
#'              list(pexp, identity, pnorm, pexp, sqrt)))
#'
#' ## Subset some variables.
#' subset(rv, c(2, 4, 3))
#' subset(rv, 5)
#' subset(rv, integer(0))
#'
#' ## This subset won't work:
#' subset(rv, c(4, 1))
#' ## But it will in a 0-truncated vine:
#' subset(trunc(rv, 0), c(4, 1))
#'
#' ## Select variables not present?
#' subset(rv, c(2, 4, 17))
#' @export
subset.rvine <- function(rv, select) {
    A <- rv$A
    copmat <- rv$copmat
    cparmat <- rv$cparmat
    p <- ncol(A)
    ntrunc <- nrow(A) - 1
    k <- length(select)
    v <- vars(rv)
    unknownvars <- setdiff(select, v)
    if (length(unknownvars) > 0)
        stop(paste("Can't subset. Variables not in vine:",
                   paste(unknownvars, collapse = ", ")))
    if (k == p) return(rv)
    if (k == 1) return(rvine(matrix(select), marg = rv$marg[[which(v == select)]]))
    if (k == 0) return(rvine(matrix(nrow=0, ncol=0)))
    ## Indices to keep (i.e. what columns of A, or which of the
    ##   ordered variables to keep?)
    ikeep <- sort(sapply(select, function(s) which(v == s)))
    ilose <- setdiff(1:p, ikeep)
    ## If it's an independence vine, it's easy to subset:
    if (ntrunc == 0) {
        A <- matrix(A[, ikeep], nrow = 1)
        blankmat <- matrix(nrow=0, ncol=k)
        return(rvine(A, marg=rv$marg[ikeep]))
    }
    ## Convert vine array to a "convenient" vine array, where the labels go on
    ##  top instead of the "skewed diagonal":
    Acon <- Atocon(A)
    ## Remove bottom rows if need be:
    if (k - 1 < ntrunc) {
        Acon <- Acon[1:k, ]
        if (!is.null(copmat)) {
            copmat <- copmat[1:(k-1), ]
            if (!is.matrix(copmat)) copmat <- matrix(copmat, nrow = 1)
        }
        if (!is.null(cparmat)) {
            cparmat <- cparmat[1:(k-1), ]
            if (!is.matrix(cparmat)) cparmat <- matrix(cparmat, nrow = 1)
        }
    }
    nrownew <- nrow(Acon)
    ## Remove columns of non-selected variables:
    Acon <- Acon[, ikeep]
    vnew <- Acon[1, ]
    if (!is.null(copmat)) copmat <- copmat[, ikeep]
    if (!is.null(cparmat)) cparmat <- cparmat[, ikeep]
    ## Go column-by-column, and gather variables that are in the selection set.
    unstrung <- integer(0)
    unstrung1 <- character(0)
    unstrung2 <- numeric(0)
    for (j in 2:k) for (i in 2:min(j, nrownew)) {
        if (Acon[i, j] %in% select) {
            unstrung <- c(unstrung, Acon[i, j])
            unstrung1 <- c(unstrung1, copmat[i-1, j])
            unstrung2 <- c(unstrung2, cparmat[i-1, j])
        }
    }
    if (length(unstrung) != nrownew * (k-1) - choose(nrownew, 2)) {
        warning(paste0("Vine can't be subsetted to selection ",
                       paste(select, collapse = ", "),
                       ". Returning NULL."))
        return(NULL)
    }
    resA <- makeuppertri(unstrung, nrow = nrownew-1, ncol = k, byRow = FALSE)
    resA <- rbind(matrix(vnew, nrow=1), resA)
    if (!is.null(copmat))
        copmat <- makeuppertri(unstrung1, nrow=nrownew-1, ncol=k, byRow=FALSE,
                               blanks = "")
    if (!is.null(cparmat)) {
        if (is.list(unstrung2)) {
            len <- sapply(unstrung2, length)
            cparmat <- makeuppertri.list(c(unstrung2, recursive=TRUE),
                                         len=len, nrow=nrownew-1, ncol=k,
                                         byRow=FALSE)
        } else {
            cparmat <- makeuppertri(unstrung2, nrow=nrownew-1, ncol=k, byRow=FALSE)
        }
    }
    ## Convert back to standard array format:
    resA <- contoA(resA)
    ## Grab the marginals:
    resmarg <- rv$marg[ikeep]
    rvine(resA, copmat, cparmat, resmarg)
}

#' @export
subset <- function(...) UseMethod("subset")


#' Subset a vine array
#'
#' Find the vine connecting a subset of variables from a bigger vine, if
#' it exists. Deprecated; use \code{\link{subset.rvine}} instead.
#'
#' @param A vine array
#' @param select Vector of variable indices in \code{diag(A)} to subset,
#' if possible. The order of the variables does not matter.
#' @details Just a technicality:
#' by saying a subset "doesn't have an existing vine", I mean that
#' a vine can't be formed using nodes and edges from
#' the original -- not that the
#' joint distribution of the selected variables can't be created from a vine
#' (so as to say, for example, that the simplifying assumption of vines
#' doesn't hold for this distribution).
#' @return Returns a vine array of the subsetted variables, with variables
#' ordered according to their order in \code{A}; or \code{NULL} if
#' the subset does not form a vine.
#' @examples
#' ## Try a D-Vine
#' rv <- rvine(CopulaModel::Dvinearray(5))
#' subset(rv, c(2, 3, 4))
#' subset(rv, c(4, 1))
#' @export
rvinesubset <- function(A, select) {
    warning("'rvinesubset' is deprecated. Please use 'subset.rvine' instead.")
    p <- ncol(A)
    ntrunc <- nrow(A) - 1
    k <- length(select)
    if (k == p) return(A)
    if (k == 1) return(matrix(select))
    if (k == 0) return(A[integer(0), integer(0)])
    vars <- varray.vars(A)
    notselect <- setdiff(vars, select)
    ## Indices to keep (i.e. what columns of A, or which of the
    ##   ordered variables to keep?)
    ikeep <- sort(sapply(select, function(s) which(vars == s)))
    ilose <- setdiff(1:p, ikeep)
    ## Convert vine array to a "convenient" vine array, where the labels go on
    ##  top instead of the "skewed diagonal":
    Acon <- A[1:ntrunc, ]
    if (!is.matrix(Acon)) Acon <- matrix(Acon, nrow = ntrunc)
    diag(Acon) <- 0
    Acon <- rbind(matrix(vars, nrow = 1), Acon)
    ## Remove bottom rows if need be:
    if (k - 1 < ntrunc) Acon <- Acon[1:k, ]
    nrownew <- nrow(Acon)
    ## Remove columns of non-selected variables:
    Acon <- Acon[, ikeep]
    ## Go column-by-column, and gather variables that are in the selection set.
    unstrung <- integer(0)
    for (j in 1:ncol(Acon)) for (i in 1:nrow(Acon)) {
        if (Acon[i, j] %in% select) unstrung <- c(unstrung, Acon[i, j])
    }
    if (length(unstrung) != nrownew * k - choose(nrownew, 2)) return(NULL)
    res <- makeuppertri(unstrung, row = nrownew, col = k,
                        incDiag = TRUE, byRow = FALSE)
    ## Convert back to standard array format:
    varsnew <- res[1, ]
    res <- res[c(2:nrownew, 1), ]
    diag(res) <- varsnew[1:nrownew]
    res[nrownew, 1:(nrownew-1)] <- 0
    res
}
