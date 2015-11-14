#' Make an upper triangular matrix
#'
#' Upper triangular matrix making made easy.
#'
#' @param entries Vector of entries for the upper triangular part of the matrix.
#' @param row Number of rows in the matrix
#' @param col Number of columns in the matrix
#' @param blanks What should go in the non-diagonal entries?
#' @param byRow Logical; read the \code{entries} by row (as if you're reading)
#' (TRUE), or read in vertically (FALSE).
#' @param incDiag Should the entries go on the diagonal too? \code{TRUE} if
#' so, \code{FALSE} if not.
#' @note Use \code{makeuppertri} to make a matrix. If you want entries to be
#' vectors (which would have to be an array with list entries), use
#' \code{makeuppertri.list}.
#' @rdname makeuppertri
#' @examples
#' ## Square matrices
#' makeuppertri(1:choose(5,2), 5, 5)
#' makeuppertri(1:choose(6,2), 5, 5, incDiag = TRUE)
#'
#' ## Not square.
#' makeuppertri(1:9, row = 2, col = 5, incDiag = TRUE)
#' makeuppertri(1:3, row = 5, col = 3)
#' @export
makeuppertri <- function(entries, row, col, blanks=0, byRow=TRUE, incDiag=FALSE){
    if (incDiag) col <- col + 1
    if (byRow) {
        comp <- matrix(blanks, row, col)
        tcomp <- t(comp)
        tcomp[t(upper.tri(comp))] <- entries
        comp <- t(tcomp)
    } else {
        comp <- matrix(blanks, row, col)
        comp[upper.tri(comp)] <- entries
    }
    if (incDiag) {
        return(comp[, 2:ncol(comp)])
    } else {
        return(comp)
    }
}

#' @param len Vector of positive integers which specify the lengths of the
#' individual vectors that are pooled in \code{entries}.
#' @examples
#' makeuppertri(1:10, 5, 5, byRow = FALSE)
#' makeuppertri.list(1:12, c(1, 10, 1), 3, 3)
#' @rdname makeuppertri
#' @export
makeuppertri.list <- function(entries, len, row, col, blanks=list(), byRow=TRUE){
    ## Put entries into list form.
    listentries <- NULL
    for (len_ in len) {
        listentries <- c(listentries, list(entries[seq(length.out = len_)]))
        ## Remove those entries that were taken (has to be done with setdiff to
        ##   account for the empty case)
        entries <- entries[setdiff(1:length(entries), seq(length.out = len_))]
    }
    ## Construct array
    if (byRow) {
        comp <- array(blanks, c(row, col))
        tcomp <- t(comp)
        tcomp[t(upper.tri(comp))] <- listentries
        comp <- t(tcomp)
    } else {
        comp <- array(blanks, c(row, col))
        comp[upper.tri(comp)] <- listentries
    }
    comp
}
