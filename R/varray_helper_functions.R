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


#' Convert Vine Array to "Convenient" Array
#'
#' A "convenient" array is achieved by moving the variables in a
#' (possibly truncated) vine array to the top row.
#'
#' @param A Vine array
#' @param Acon A convenient vine array
#' @details
#' \code{Atocon} converts a vine array to a convenient vine array.
#' \code{contoA} converts a convenient vine array to a vine array.
#'
#' These converters are intended to be used internally.
#' @note This form of array is "convenient" because truncated vines can
#' be respresented by truncating the bottom rows of the matrix.
#' @return A vine array or convenient vine array
#' @rdname A_Acon_convert
#' @export
Atocon <- function(A) {
    ntrunc <- nrow(A) - 1
    if (ntrunc == 0) return(A)
    d <- ncol(A)
    v <- c(diag(A), A[ntrunc + 1, ntrunc+1+seq_len(d-ntrunc-1)])
    Acon <- A[1:ntrunc, ]
    if (!is.matrix(Acon)) Acon <- matrix(Acon, nrow = ntrunc)
    diag(Acon) <- 0
    rbind(matrix(v, nrow = 1), Acon)
}

#' @rdname A_Acon_convert
#' @export
contoA <- function(Acon) {
    ntrunc <- nrow(Acon) - 1
    if (ntrunc == 0) return(Acon)
    d <- ncol(Acon)
    v <- Acon[1, ]
    A <- Acon[-1, ]
    if (!is.matrix(A)) A <- matrix(A, nrow = 1)
    A <- rbind(A, matrix(c(rep(0, ntrunc), v[(ntrunc+1):d]), nrow = 1))
    diag(A) <- v[1:(ntrunc+1)]
    A
}
