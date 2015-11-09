#' Subset a vine array
#'
#' Find the vine connecting a subset of variables from a bigger vine, if
#' it exists.
#'
#' @param A vine array
#' @param select Vector of variable indices in \code{diag(A)} to subset,
#' if possible.
#' @details Just a technicality:
#' by saying a subset "doesn't have an existing vine", I mean that
#' a vine can't be formed using nodes and edges from
#' the original -- not that the
#' joint distribution of the selected variables can't be created from a vine
#' (so as to say, for example, that the simplifying assumption of vines
#' doesn't hold for this distribution).
#' @return Returns a vine array of the subsetted variables, or \code{NULL} if
#' the subset don't form a vine.
#' @examples
#' ## Try a D-Vine
#' A <- CopulaModel::Dvinearray(5)
#' rvinesubset(A, c(2, 3, 4))
#' rvinesubset(A, c(4, 1))
#' ## The numbering doesn't have to be 1,2,3,...
#' A <- makeuppertri(2 + c(1,1,1,1,4, 2,2,3,1, 3,2,3, 4,2, 5), 6, 6)[1:5, 2:6]
#' rvinesubset(A, 2 + c(1,3,4))
#' @export
rvinesubset <- function(A, select) {
    p <- ncol(A)
    k <- length(select)
    if (k == p) return(A)
    if (k == 1) return(matrix(select))
    diagA <- diag(A)
    notselect <- setdiff(diag(A), select)
    ## Indices to keep (i.e. what columns of A, or which of the
    ##   ordered variables to keep?)
    ikeep <- sapply(select, function(s) which(diagA == s))
    ilose <- setdiff(1:p, ikeep)
    ## Select entries to remove from the vine array (as TRUE entries).
    #### Bottom-left zeroes need removal.
    removal <- lower.tri(diag(p))
    #### "Extraneous" area needs removal (i.e. higher-level trees
    ####   that are impossible with this selection)
    removal[k:p, k:p] <- removal[k:p, k:p] | upper.tri(diag(p-k+1))
    #### Columns where our selection is not on the diagonal need removing:
    removal[, ilose] <- TRUE
    #### Select remaining variables that need removal.
    removal[, ikeep] <- removal[, ikeep] |
        apply(A[, ikeep], 1:2, function(t) t %in% notselect)
    ## There should be (k+1) choose 2 variable indices remaining if the
    ##  subsetted vine "exists".
    remain <- A[!removal]  # A vector.
    if (length(remain) == choose(k+1, 2)) {
        subA <- makeuppertri(remain, k+1, k+1, byRow = FALSE)
        subA <- subA[, -1]
        subA <- subA[-(k+1), ]
    } else {
        subA <- NULL
    }
    subA
}
