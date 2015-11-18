#' "Re-leaf" a Vine Array
#'
#' Convert a vine array so that a variable(s) of your choice appears last
#' (bottom-right corner) of the vine array -- if possible.
#'
#' @param A Vine array matrix, possibly truncated.
#' @param leaf Vector of variables that you want to
#' get new vine array for, if such a vine array exists.
#' @return A vine array (matrix) with \code{leaf} at the end, if such a
#' vine array exists; \code{NULL} otherwise. If \code{length(leaf) > 1},
#' a list of such output is returned, corresponding to the entries
#' in \code{leaf}.
#' @examples
#' A <- makeuppertri(c(4,4,4,4,6,5,
#'                     6,6,6,4,4,
#'                     1,1,1,6,
#'                     5,5,1,
#'                     2,2,
#'                     3), 6, 6, incDiag=T)
#' releaf.varray(truncvarray(A, 2))
#' releaf.varray(truncvarray(A, 3), leaf = 3:6)
#' releaf.varray(A, leaf = 5)
#'
#' ## The function doesn't require 1:ncol(A) as variables.
#' A <- truncvarray(CopulaModel::Dvinearray(5), 2)
#' A <- apply(A, 1:2, function(i) if (i==0) "" else letters[i])
#' releaf.varray(A)
#' @export
releaf.varray <- function(A, leaf=varray.vars(A)) {
    if (!is.matrix(A)) return(NULL)
    if (length(leaf) == 0) return(list())
    if (ncol(A) == 0) return(NULL) # Only after checking length(leaf) > 0.
    d <- ncol(A)
    ntrunc <- nrow(A) - 1
    varsset <- varray.vars(A)
    ## Case when d=2
    if (d == 2) {
        res <- lapply(leaf, function(l) {
            if (A[2,2] == l) {
                A
            } else {
                matrix(c(l, 0, l, A[2,2]), ncol = 2)
            }
        })
        if (length(res) == 1) res <- res[[1]]
        return(res)
    }
    ## It's easiest to work with arrays truncated < d-1. So make d-2 the max.
    ##  After all, it's easy to go from d-2 back to the full d-1, which
    ##  we'll do later.
    putback <- FALSE
    if (ntrunc == d-1) {
        putback <- TRUE
        ntrunc <- d-2
        A <- truncvarray(A, ntrunc)
    }
    ## Centralize the vine array:
    Ac <- center.varray(A)
    leftAc <- Ac[, 1:ntrunc]
    if (!is.matrix(leftAc)) leftAc <- matrix(leftAc, ncol = ntrunc)
    ## Identify truncated region to the right of column ntrunc. They're
    ##  candidate columns/variables.
    vars <- varray.vars(Ac)
    candvars <- vars[(ntrunc+1):d]
    candbranches <- Ac[1:ntrunc, (ntrunc+1):d]
    candcols <- Ac[1:(ntrunc+1), (ntrunc+1):d]
    if (!is.matrix(candbranches))
        candbranches <- matrix(candbranches, ncol = ntrunc)
    ncols <- ncol(candcols)
    ## Variables that won't work as leaves:
    bad <- unique(as.vector(candbranches))
    works <- !(leaf %in% bad)
    ## Re-order the vine arrays:
    res <- list()
    for (i in 1:length(leaf)) {
        if (works[i]) {
            ## Leaf is valid.
            col <- which(candvars == leaf[i])
            recol <- c((1:ncols)[-col], col)
            res[[i]] <- cbind(leftAc, candcols[, recol])
            if (putback) {
                res[[i]] <- rbind(res[[i]], matrix(0, ncol = d))
                res[[i]][d,d] <- res[[i]][d-1, d]
                res[[i]][d-1,d] <- setdiff(varsset, res[[i]][, d])
            }
        } else {
            ## Leaf is invalid.
            res <- c(res, list(NULL))
        }
    }
    names(res) <- leaf
    if (length(res) == 1) res <- res[[1]]
    res
}
