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

