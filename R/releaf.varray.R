#' "Re-leaf" a Vine Array
#'
#' Convert a vine array so that a variable(s) of your choice appears last
#' (bottom-right corner) of the vine array -- if possible.
#'
#' @param A Vine array matrix, possibly truncated. Must be in natural order --
#' if truncated, that means \code{A[, 1+0:nrow(A)]} must be in natural order.
#' @param leaf Vector of variables that you want to
#' get new vine array for, if such a vine array exists.
#' @return A vine array (matrix) with \code{leaf} at the end, if such a
#' vine array exists; \code{NULL} otherwise. If \code{length(leaf) > 1},
#' a list of such output is returned, corresponding to the entries
#' in \code{leaf}.
#' @note
#' There's no functionality yet to convert a truncated vine so that it's leftern
#' part is in natural order. But if you truncate a vine in natural order, it'll
#' be fine. Use \code{\link{CopulaModel::varray2NO}} for converting a complete
#' vine array to natural order.
#' @examples
#' A <- makeuppertri(c(4,4,4,4,6,5,
#'                     6,6,6,4,4,
#'                     1,1,1,6,
#'                     5,5,1,
#'                     2,2,
#'                     3), 6, 6, incDiag=T)
#' releaf.varray(trunc.varray(A, 2))
#' releaf.varray(trunc.varray(A, 3), leaf = 3:6)
#' releaf.varray(A, leaf = 5)
#' @export
releaf.varray <- function(A, leaf=diag(A)) {
    if (!is.matrix(A)) return(NULL)
    if (length(leaf) == 0) return(logical(0))
    if (ncol(A) == 0) return(NULL) # Only after checking length(leaf) > 0.
    d <- ncol(A)
    ntrunc <- nrow(A) - 1
    # if (ntrunc == d-1) ntrunc <- d-2
    ## Convert to natural order:
    ANO <- CopulaModel::varray2NO(A)$NOa
    diagNO <- diag(ANO)
    leftANO <- ANO[1:(ntrunc+1), 1:ntrunc]
    if (!is.matrix(leftANO)) leftANO <- matrix(leftANO, ncol = ntrunc)
    ## Identify truncated region to the right of column ntrunc. They're
    ##  candidate columns/variables.
    candvars <- diagNO[(ntrunc+1):d]
    candbranches <- ANO[1:ntrunc, (ntrunc+1):d]
    if (!is.matrix(candbranches)) candbranches <- matrix(candbranches, ncol = ntrunc)
    candcols <- rbind(candbranches, matrix(candvars, nrow = 1))
    ncols <- ncol(candcols)
    ## Variables that won't work as leaves:
    bad <- unique(as.vector(candbranches))
    yesorno <- !(leaf %in% bad)
    ## Re-order the vine arrays:
    res <- list()
    for (i in 1:length(leaf)) {
        if (yesorno[i]) {
            ## Leaf is valid.
            col <- which(candvars == leaf[i])
            recol <- c((1:ncols)[-col], col)
            res[[i]] <- cbind(leftANO, candcols[, recol])
            if (ntrunc == d-2) {
                res[[i]] <- rbind(res[[i]], matrix(0, ncol = d))
                res[[i]][d,d] <- res[[i]][d-1, d]
                res[[i]][d-1,d] <- setdiff(diagNO, res[[i]][, d])
            }
        } else {
            ## Leaf is invalid.
            res[[i]] <- NULL
        }
    }
    names(res) <- leaf
    if (length(res) == 1) res <- res[[1]]
    res
}
