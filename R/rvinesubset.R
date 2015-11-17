#' Subset a vine array
#'
#' Find the vine connecting a subset of variables from a bigger vine, if
#' it exists.
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
#' A <- CopulaModel::Dvinearray(5)
#' rvinesubset(A, c(2, 3, 4))
#' rvinesubset(A, c(4, 1))
#' @export
rvinesubset <- function(A, select) {
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
