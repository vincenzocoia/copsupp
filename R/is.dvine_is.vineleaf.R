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

#' Are variables Vine Leaves?
#'
#' Check whether some variables on a vine are "leaves" -- that is, only
#' connect to one other variable on the first tree.
#'
#' @param A Vine array matrix.
#' @param select Vector of vine variables that you want to query.
#' @return Vector of length \code{length(select)} of logicals.
#' @examples
#' A <- CopulaModel::Dvinearray(5)
#' is.vineleaf(A, 5)
#' A <- rvinesubset(A, 2:4)
#' is.vineleaf(A, 3:4)
#' is.vineleaf(A, integer(0))
#' @export
is.vineleaf <- function(A, select=diag(A)) {
    if (!is.matrix(A)) return(NULL)
    if (length(select) == 0) return(logical(0))
    if (ncol(A) == 0) return(NULL) # Only after checking length(select) > 0.
    ## Get nodes on tree 1:
    nodes1 <- A[1, -1]
    nodes2 <- diag(A)[-1]
    freq <- table(c(nodes1, nodes2))
    freq.select <- freq[as.character(select)]
    res <- freq.select == 1
    attributes(res) <- NULL
    res
}
