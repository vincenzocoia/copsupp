#' Is a Vine Array a D-Vine?
#'
#' @param A Vine array matrix.
#' @return Logical -- \code{TRUE} if \code{A} is a D-vine. \code{FALSE} if not.
#' If \code{A} is not a matrix, or has \code{integer(0)} columns,
#' it'll return \code{NULL}.
#' @examples
#' is.dvine(CopulaModel::Dvinearray(5))
#' is.dvine(CopulaModel::Cvinearray(6))
#' is.dvine("hello")
#' @export
is.dvine <- function(A) {
    if (!is.matrix(A)) return(NULL)
    p <- ncol(A)
    if (p == 1) return(TRUE)
    if (p == 0) return(NULL)
    ## In tree 1, nodes should appear maximum of two times.
    nodes1 <- A[1, -1]
    nodes2 <- diag(A)[-1]
    max(table(c(nodes1, nodes2))) <= 2
}

#' "Re-leaf" a Vine Array
#'
#' Convert a vine array so that a variable(s) of your choice appears last
#' (bottom-right corner) of the vine array -- if possible.
#'
#' @param A Vine array matrix.
#' @param ntrunc Truncation integer, between 1 and \code{ncol(A)-1}.
#' @param leaves Vector of variables in \code{diag(A)} that you want to
#' get new vine array for, if such a vine array exists.
#' @return A vine array (matrix) with \code{leaves} at the end, if such a
#' vine array exists; \code{NULL} otherwise. If \code{length(leaves) > 1},
#' a list of such output is returned, corresponding to the entries
#' in \code{leaves}
#' @note
#' When the truncation is less than \code{ncol(A) - 2}, the array matrices
#' that are returned are not \code{ncol(A) x ncol(A)} -- they're
#' \code{(ntrunc+1) x ncol(A)}, achieved by collapsing the variables on
#' the diagonal upwards.
#' @examples
#' A <- makeuppertri(c(1,1,1,1,2,4,
#'                     2,2,2,1,1,
#'                     3,3,3,2,
#'                     4,4,3,
#'                     5,5,
#'                     6), 7, 7)[1:6, 2:7]
#' releafvarray(A, ntrunc = 2)
#' releafvarray(A, ntrunc = 3, leaves = 3:6)
#' releafvarray(A, leaves = 5)
#' @export
releafvarray <- function(A, ntrunc = ncol(A)-1, leaves=diag(A)) {
    if (!is.matrix(A)) return(NULL)
    if (length(leaves) == 0) return(logical(0))
    if (ncol(A) == 0) return(NULL) # Only after checking length(leaves) > 0.
    d <- ncol(A)
    if (ntrunc == d-1) ntrunc <- d-2
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
    yesorno <- !(leaves %in% bad)
    ## Re-order the vine arrays:
    res <- list()
    for (i in 1:length(leaves)) {
        if (yesorno[i]) {
            ## Leaf is valid.
            col <- which(candvars == leaves[i])
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
    names(res) <- leaves
    if (length(res) == 1) res <- res[[1]]
    res
}
