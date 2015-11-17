#' Rearrange the Matrices following Rearrangement of Vine Array
#'
#' If you changed the order of your vine array, this function will
#' convert other matrices associated with the original vine array
#' (such as copula matrix or copula parameter matrix)
#' into new matrices associated with the new vine array.
#'
#' @param mat Old matrix to rearrange, such as a copula matrix or copula
#' parameter matrix. Should be upper-triangular, with nrows = truncation number,
#' ncols = vine dimension.
#' @param Anew The new vine array index, possibly truncated.
#' @param Aold The old vine array index, possibly truncated.
#' @examples
#' ## Originally...
#' A <- truncvarray(Cvinearray(4), 2)
#' copmat <- makeuppertri(c("gum", "gal", "bvtcop",
#'                          "bvncop", "frk"), row = 2, col = 4, blanks = "")
#' cparmat <- makeuppertri.list(c(2, 3, 0.9, 4, 0.1, 0.5),
#'                              len = c(1,1,2,1,1), row = 2, col = 4)
#'
#' ## Obtain new vine array by some means:
#' Anew <- center.varray(A)
#' ## Get new matrices:
#' reform.copmat(copmat, Anew, A)
#' reform.copmat(cparmat, Anew, A)
#' @export
reform.copmat <- function(mat, Anew, Aold) {
    ntrunc <- nrow(Anew) - 1
    d <- ncol(Anew)
    ## Convert to convenient format:
    Anew <- Atocon(Anew)
    Aold <- Atocon(Aold)
    ## Start new matrix:
    newmat <- matrix(mat[nrow(mat), 1], nrow = nrow(mat), ncol = ncol(mat))
    for (i in 1+1:ntrunc) for (j in i:d) {
        ## Fill in newmat[i-1, j]. First find the variable and conditioning pair
        pair <- c(Anew[1, j], Anew[i, j])
        cond <- Anew[1+seq_len(i-2), j]
        ## Which entry was that in the old A? Collect the column number.
        for (jold in i:d) {
            pairold <- c(Aold[1, jold], Aold[i, jold])
            if (i == 2) {
                if (all(sort(pair) == sort(pairold))) break
            } else {
                condold <- Aold[2:(i-1), jold]
                if (all(c(sort(cond) == sort(condold),
                          sort(pair) == sort(pairold)))) break
            }
        }
        ## Grab the corresponding entry in mat:
        newmat[i-1, j] <- mat[i-1, jold]
    }
    newmat
}
