#' Relabel Variables in a Vine Array
#'
#' @param A Vine array. Could be truncated.
#' @param labs Vector of new labels. The order of the labels correspond to
#' the order of the variables in the vine array \code{A}.
#' @details Identical to \code{CopulaModel::varrayperm} but allows for
#' the posibility that \code{A} is not square.
#' @return A relabelled vine array matrix.
#' @examples
#' (A <- trunc.varray(CopulaModel::Cvinearray(5), 2))
#' relabel.varray(A, c(3, 2, 1, 5, 4))
#' @export
relabel.varray <- function(A, labs = 1:ncol(A)) {
    p <- ncol(A)
    r <- nrow(A)
    labs_orig <- varray.vars(A)
    wch_lab <- invert.perm(labs_orig)
    for (row in 1:r) {
        for (col in row:p) {
            A[row, col] <- labs[wch_lab[A[row, col]]]
        }
    }
    A
}

#' Invert a permutation
#'
#' For a permutation of a set of integers \code{1, 2, ...,p}, finds the
#' inverse permutation using the \code{\link{which}} function.
#'
#' @param perm Vector of integers in \code{{1:length(perm)}}
#' @return A vector of length \code{length(perm)} of the inverse permutation.
#' @note This function won't check whether the integers you input are
#' in the set \code{{1:length(perm)}}, but allows for length-0 entry.
#' @examples
#' perm <- c(5, 1, 2, 3, 4)
#' (perminv <- invert.perm(perm))
#' perm[perminv]
#' perminv[perm]
#'
#' ## The zero case:
#' invert.perm(integer(0))
#' @export
invert.perm <- function(perm) {
    p <- length(perm)
    if (p <= 1) return(perm)
    sapply(1:p, function(i) which(perm == i))
}

#' Truncate a Vine Array
#'
#' Truncates a vine array, collapsing the variables upwards.
#'
#' @param A A vine array, possibly truncated.
#' @param ntrunc Integer; truncation level
#' @details If \code{ntrunc >= nrow(A) - 1}, the original vine array is
#' returned. Otherwise, a truncated vine array with \code{ntrunc + 1} rows is
#' returned. The variables are listed along the initial diagonal of the vine
#' array, then continue along the bottom row.
#' @examples
#' (A <- CopulaModel::Dvinearray(6))
#' (A <- trunc.varray(A, 3))
#' (A <- relabel.varray(A, c(6, 2, 4, 3, 1, 5)))
#' trunc.varray(A, 2)
#' @export
trunc.varray <- function(A, ntrunc) {
    p <- ncol(A)
    r <- nrow(A)
    if (ntrunc > r-2) return(A)
    vars <- varray.vars(A)
    A <- A[1:ntrunc, 1:p]
    vars[1:ntrunc] <- 0
    rbind(A, matrix(vars, nrow=1))
}

#' Extract Variables in a Vine Array
#'
#' Extract variables in a vine array, possibly truncated, in the order
#' of the vine array.
#'
#' @param A Vine array, possibly truncated.
#' @return Vector of vine variables.
#' @examples
#' A <- CopulaModel::Dvinearray(5)
#' A <- relabel.varray(A, c(5, 2, 4, 3, 1))
#' varray.vars(A)
#'
#' A <- trunc.varray(A, 2)
#' varray.vars(A)
#' @export
varray.vars <- function(A) {
    p <- ncol(A)
    r <- nrow(A)
    firstvars <- diag(A)
    secondvars <- A[r, r+seq_len(p-r)]
    c(firstvars, secondvars)
}

#' Center a Vine Array
#'
#' Converts a vine array \code{A} so that the first variables (up to
#' truncation) are not leaves. So, a slightly weaker condition than
#' natural order.
#'
#' @param A Vine array, possibly truncated.
#' @details For a \code{t}-truncated vine array \code{(t < ncol(A)-1)},
#' the vine array is re-ordered so that the first \code{t} variables
#' introduced in the outputted array are not leaves.
#'
#' If \code{t = ncol(A)-1}, then the entered vine isn't truncated, and the
#' first \code{t-1} variables in the outputted array are not leaves (in fact,
#' a natural order array is outputted, since it satisfies that requirement).
#' @export
center.varray <- function(A) {
    ntrunc <- nrow(A) - 1
    d <- ncol(A)
    if (ntrunc == d-1) return(CopulaModel::varray2NO(A)$NOa)
    ## Initiate the final array as (ntrunc+1)x(ntrunc+1) array using variables
    ##  in A[, d], with A[d,d] going at the end.
}

