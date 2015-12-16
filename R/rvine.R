#' Create an object of class `rvine`
#'
#' @param A Vine array (matrix), possibly truncated.
#' @param copmat Upper triangular \code{(nrow(A)-1) x ncol(A)} matrix of copula model names.
#' @param cparmat Upper triangular \code{(nrow(A)-1) x ncol(A)} matrix of
#' parameters taken by the copulas in \code{copmat}.
#' @param marg List of \code{ncol(A)} vectorized marginal distribution functions
#' with entries corresponding to the columns in \code{A}.
#' @return An object of class `rvine`, which is a named list of the arguments.
#' @examples
#' A <- CopulaModel::Dvinearray(4)
#' copmat <- makeuppertri(c("gum", "bvtcop", "mtcj",
#'                          "frk", "indepcop",
#'                          "frk"), 3, 4, blanks = "")
#' cparmat <- makeuppertri.list(c(2, 0.5, 4, 2,
#'                                1,
#'                                1), len = c(1,2,1,1,0,1), nrow = 3, ncol = 4)
#' rvine(A, copmat, cparmat)
#' @export
rvine <- function(A, copmat = NULL, cparmat = NULL, marg = NULL) {
    res <- list(A = A,
                copmat = copmat,
                cparmat = cparmat,
                marg = marg)
    class(res) <- "rvine"
    res
}
