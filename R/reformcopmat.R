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
#' @param Anew The new vine array index, possibly truncated. Could have
#' less variables and more truncation than \code{Aold}.
#' @param Aold The old vine array index, possibly truncated.
#' @note Make sure that the variables in \code{Aold}
#' and \code{Anew} have the same labels.
#' @return
#' A matrix that's possibly subsetted and/or rearranged from \code{mat}.
#' @examples
#' ## Originally...
#' G <- AtoG(Cvinearray(4))[1:3, ]
#' copmat <- makeuppertri(c("gum", "gal", "bvtcop",
#'                          "bvncop", "frk"), row = 2, col = 4, blanks = "")
#' cparmat <- makeuppertri.list(c(2, 3, 0.9, 4, 0.1, 0.5),
#'                              len = c(1,1,2,1,1), row = 2, col = 4)
#'
#' ## Obtain new vine array by some means:
#' Gnew <- centervarray(G)
#' ## Get new matrices:
#' reformcopmat(copmat, Gnew, G)
#' reformcopmat(cparmat, Gnew, G)
#'
#' ## Try changing the dimension of G
#' Gnew <- subsetvarray(G[1:2, ], 1:3)
#' reformcopmat(copmat, Gnew, G)
#' reformcopmat(cparmat, Gnew, G)
#' @export
reformcopmat <- function(mat, Gnew, Gold) {
    ntrunc <- nrow(Gnew) - 1
    d <- ncol(Gnew)
    if (ntrunc == 0) return(matrix(nrow=0, ncol=d))
    dold <- ncol(Gold)
    ## Start new matrix:
    entry <- mat[nrow(mat), 1]  # "" or 0?
    newmat <- matrix(entry, nrow = ntrunc, ncol = d)
    for (i in 1+1:ntrunc) for (j in i:d) {
        ## Fill in newmat[i-1, j]. First find the variable and conditioning pair
        pair <- c(Gnew[1, j], Gnew[i, j])
        cond <- Gnew[1+seq_len(i-2), j]
        ## Which entry was that in the old A? Collect the column number.
        for (jold in i:dold) {
            pairold <- c(Gold[1, jold], Gold[i, jold])
            if (i == 2) {
                if (all(sort(pair) == sort(pairold))) break
            } else {
                condold <- Gold[2:(i-1), jold]
                if (all(c(sort(cond) == sort(condold),
                          sort(pair) == sort(pairold)))) break
            }
        }
        ## Grab the corresponding entry in mat:
        newmat[i-1, j] <- mat[i-1, jold]
    }
    newmat
}
