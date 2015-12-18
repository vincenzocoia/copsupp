#' Relabel Variables in a Regular Vine
#'
#' @param rv Object of type 'rvine'.
#' @param labs Vector of new labels. The order of the labels correspond to
#' the order of the variables in the vine array \code{A}.
#' @details Similar to \code{CopulaModel::varrayperm} but allows for
#' non-square vine arrays, as well as labels outside
#' of the set \code{{1:ncol(rv$A)}}.
#' @return An object of type 'rvine' with a relabelled vine array matrix.
#' @examples
#' (rv <- rvine(truncvarray(CopulaModel::Cvinearray(5), 2)))
#' relabel(rv, c(3, 2, 1, 5, 4))
#'
#' ## Labels don't need to be in the set {1:ncol(A)}.
#' A <- CopulaModel::Dvinearray(4)
#' A <- apply(A, 1:2, function(i) if (i==0) "" else letters[i])
#' relabel.varray(rvine(A), c("this", "is", "a", "demo"))
#' relabel.varray(rvine(A))
#' @export
relabel.rvine <- function(rv, labs = 1:ncol(rv$A)) {
    A <- rv$A
    d <- ncol(A)
    r <- nrow(A)
    labs_orig <- vars(rv)
    isnumeric <- is.numeric(labs)
    ## Map original label to order
    map2order <- function(labo) which(labs_orig == labo)
    ## Map original labels to new labels
    for (row in 1:r) {
        for (col in row:d) {
            A[row, col] <- labs[map2order(A[row, col])]
        }
    }
    ## This is a little extensive, to allow for non-integers...
    ##   but why not include it?
    if (isnumeric){
        A[lower.tri(A)] <- 0
        A <- apply(A, 1:2, as.numeric)
    }
    rvine(A, rv$copmat, rv$cparmat, rv$marg)
}

#' @export
relabel <- function(...) UseMethod("relabel")

#' Relabel Variables in a Vine Array
#'
#' Deprecated; use \code{\link{relabel.rvine}} instead.
#'
#' @param A Vine array. Could be truncated.
#' @param labs Vector of new labels. The order of the labels correspond to
#' the order of the variables in the vine array \code{A}.
#' @details Similar to \code{CopulaModel::varrayperm} but allows for
#' the posibility that \code{A} is not square, as well as labels outside
#' of the set \code{{1:ncol(A)}}.
#' @return A relabelled vine array matrix.
#' @examples
#' (A <- truncvarray(CopulaModel::Cvinearray(5), 2))
#' relabel.varray(A, c(3, 2, 1, 5, 4))
#'
#' ## Labels don't need to be in the set {1:ncol(A)}.
#' A <- CopulaModel::Dvinearray(4)
#' A <- apply(A, 1:2, function(i) if (i==0) "" else letters[i])
#' relabel.varray(A, c("this", "is", "a", "demo"))
#' relabel.varray(A)
#' @export
relabel.varray <- function(A, labs = 1:ncol(A)) {
    warning("Function 'relabel.varray' is deprecated. Please use 'relabel.rvine'.")
    d <- ncol(A)
    r <- nrow(A)
    labs_orig <- varray.vars(A)
    isnumeric <- is.numeric(labs)
    ## Map original label to order
    map2order <- function(labo) which(labs_orig == labo)
    ## Map original labels to new labels
    for (row in 1:r) {
        for (col in row:d) {
            A[row, col] <- labs[map2order(A[row, col])]
        }
    }
    ## This is a little extensive, to allow for non-integers...
    ##   but why not include it?
    if (isnumeric){
        A[lower.tri(A)] <- 0
        return(apply(A, 1:2, as.numeric))
    } else {
        return(A)
    }
}
