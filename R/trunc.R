#' Truncate a Regular Vine
#'
#' Truncates a regular vine.
#'
#' @param rv A regular vine object.
#' @param ntrunc Integer; truncation level.
#' @return If \code{ntrunc >= nrow(A) - 1}, the original vine is
#' returned. Otherwise, an 'rvine' object is returned
#' with the truncation implemented.
#' @note The variables are listed along the initial diagonal of the vine
#' array, then continue along the bottom row.
#' @examples
#' rv <- rvine(CopulaModel::Dvinearray(6))
#' trunc(rv, 3)
#' trunc(rv, 100)
#' trunc(rv, 0)
#' @export
trunc.rvine <- function(rv, ntrunc) {
    A <- rv$A
    copmat <- rv$copmat
    cparmat <- rv$cparmat
    d <- ncol(A)
    r <- nrow(A)
    if (ntrunc > r-2) return(rv)
    if (ntrunc == 0) return(rvine(matrix(vars(rv), nrow=1), marg=rv$marg))
    Acon <- Atocon(A)
    Acon <- Acon[1:(ntrunc+1), ]
    A <- contoA(Acon)
    if (!is.matrix(Acon)) Acon <- matrix(Acon, nrow = 1)
    if (!is.null(copmat)) {
        copmat <- copmat[seq_len(ntrunc), ]
        if (!is.matrix(copmat)) copmat <- matrix(copmat, nrow = 1)
    }
    if (!is.null(cparmat)) {
        cparmat <- cparmat[seq_len(ntrunc), ]
        if (!is.matrix(cparmat)) cparmat <- matrix(cparmat, nrow = 1)
    }
    rvine(A, copmat, cparmat, rv$marg)
}

#' @export
trunc <- function(...) UseMethod("trunc")

#' Truncate a Vine Array
#'
#' Truncates a vine array, collapsing the variables upwards. Deprecated;
#' please use \code{\link{trunc.rvine}} instead.
#'
#' @param A A vine array, possibly truncated.
#' @param ntrunc Integer; truncation level
#' @return If \code{ntrunc >= nrow(A) - 1}, the original vine array is
#' returned. Otherwise, a truncated vine array with \code{ntrunc + 1} rows is
#' returned.
#' @note The variables are listed along the initial diagonal of the vine
#' array, then continue along the bottom row.
#' @examples
#' (A <- CopulaModel::Dvinearray(6))
#' (A <- truncvarray(A, 3))
#' (A <- relabel.varray(A, c(6, 2, 4, 3, 1, 5)))
#' truncvarray(A, 2)
#' @export
truncvarray <- function(A, ntrunc) {
    d <- ncol(A)
    r <- nrow(A)
    if (ntrunc > r-2) return(A)
    vars <- varray.vars(A)
    A <- A[1:ntrunc, 1:d]
    if (!is.matrix(A)) A <- matrix(A, ncol = d)
    vars[1:ntrunc] <- 0
    rbind(A, matrix(vars, nrow=1))
}
