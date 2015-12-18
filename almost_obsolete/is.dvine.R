#' Is a Vine Array a D-Vine?
#'
#' @param A Vine array, possibly truncated.
#' @return Logical -- \code{TRUE} if \code{A} is a D-vine. \code{FALSE} if not.
#' If \code{A} is not a matrix, or has \code{integer(0)} columns,
#' it'll return \code{NULL}.
#' @examples
#' A <- CopulaModel::Dvinearray(5)
#' is.dvine(A)
#' A <- truncvarray(CopulaModel::Cvinearray(6), 2)
#' is.dvine(A)
#' is.dvine("hello")
#' @export
is.dvine <- function(A) {
    if (!is.matrix(A)) return(NULL)
    p <- ncol(A)
    if (p == 1) return(TRUE)
    if (p == 0) return(NULL)
    ## In tree 1, nodes should appear maximum of two times.
    nodes1 <- varray.vars(A)[-1]
    nodes2 <- A[1, -1]
    max(table(c(nodes1, nodes2))) <= 2
}

