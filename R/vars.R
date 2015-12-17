#' Extract Variables from a Regular Vine
#'
#' Extract variables from a regular vine, possibly truncated, in the order
#' of the vine array.
#'
#' @param rv A regular vine object.
#' @return Vector of vine variables.
#' @examples
#' rv <- rvine(CopulaModel::Dvinearray(5))
#' rv <- relabel(A, c(5, 2, 4, 3, 1))
#' vars(rv)
#'
#' rv <- trunc(rv, 2)
#' vars(rv)
#' @export
vars.rvine <- function(rv) {
    A <- rv$A
    p <- ncol(A)
    if (p == 0) return(integer(0))
    r <- nrow(A)
    firstvars <- diag(A)
    secondvars <- A[r, r+seq_len(p-r)]
    c(firstvars, secondvars)
}

#' @export
vars <- function(...) UseMethod("vars")

#' Extract Variables in a Vine Array
#'
#' Extract variables in a vine array, possibly truncated, in the order
#' of the vine array. Deprecated; use \code{\link{vars.rvine}} instead.
#'
#' @param A Vine array, possibly truncated.
#' @return Vector of vine variables.
#' @examples
#' A <- CopulaModel::Dvinearray(5)
#' A <- relabel.varray(A, c(5, 2, 4, 3, 1))
#' varray.vars(A)
#'
#' A <- truncvarray(A, 2)
#' varray.vars(A)
#' @export
varray.vars <- function(A) {
    warning("'varray.vars' is deprecated. Please use 'vars.rvine'.")
    p <- ncol(A)
    r <- nrow(A)
    firstvars <- diag(A)
    secondvars <- A[r, r+seq_len(p-r)]
    c(firstvars, secondvars)
}
