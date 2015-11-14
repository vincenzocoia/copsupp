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
#' @note Use \code{makeuppertri} to make a matrix. If you want entries to be
#' vectors (which would have to be an array with list entries), use
#' \code{makeuppertri.list}.
#' @rdname makeuppertri
#' @export
makeuppertri <- function(entries, row, col, blanks=0, byRow=TRUE){
    if (byRow) {
        comp <- matrix(blanks, row, col)
        tcomp <- t(comp)
        tcomp[t(upper.tri(comp))] <- entries
        comp <- t(tcomp)
    } else {
        comp <- matrix(blanks, row, col)
        comp[upper.tri(comp)] <- entries
    }
    comp
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
