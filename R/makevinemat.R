#' Construct a Matrix for a Vine
#'
#' Use this function to construct a vine array, copula matrix, or copula parameter
#' matrix. Do so by specifying layers in the vine (columns in the matrix)
#' or trees in the vine (rows in the matrix). A full matrix need not be
#' constructed.
#'
#' @param ... Vectors or lists to construct the matrix with.
#' @param bylayer \code{TRUE} if \code{...} are layers (columns in the outputted
#' matrix). \code{FALSE} if they are trees (rows in the outputted matrix).
#' @param zerocol Should there be an "empty" column appended to the left
#' of the outputted matrix? Useful for making copula matrices, for example.
#' @return If \code{bylayer} is \code{TRUE}, the returned matrix places
#' entries in "..." as
#' columns from left to right, as if building a vine by adding layers.
#' Empty portions of the columns are put at the
#' bottom of the column.
#'
#' If \code{bylayer} is \code{FALSE}, the returned matrix places
#' entries in "..." as
#' rows from top to bottom, as if building a vine by specifying tree depths.
#' Empty portions of the rows are put to the
#' left of the row.
#' @details
#' If you want an entry of the matrix to contain something other than a vector
#' of length 1 (for example, copula parameters of dimension different from 1),
#' specify the entries within a list in \code{...}.
#'
#' "Empty" entries are either \code{list(NULL)} if the matrix entries aren't
#' all vectors of length 1, \code{0} if the matrix is numeric, or \code{""}
#' if the matrix is of type "character".
#' @note If you intend to build a vine array with this function, you'll get
#' "tree 0" (the variable nodes) in row 1, instead of the diagonal. To convert
#' to the "traditional" form, use \code{\link{contoA}}.
#' @examples
#' ## Vine array
#' makevinemat(1, c(2,1), c(3, 1, 2), c(4, 1, 2, 3))
#' makevinemat(1, c(2,1), c(3, 1, 2), c(4, 1), c(5, 1:4))
#'
#' ## Two layers of a vine array
#' makevinemat(c(4,1), c(5, 1:4))
#'
#' ## Copula Matrix
#' makevinemat(c("gum", "gum", "bvtcop"), c("frk", "", "joe"),
#'             bylayer=F, zerocol=T)
#'
#' ## Copula Parameter Matrix
#' makevinemat(list(2.4, 3.8, c(0.6, 7)), c(4.5, 2.3), bylayer=F, zerocol=T)
#' makevinemat(1:3, 4:5, bylayer=F, zerocol=T)
#' @export
makevinemat <- function(..., bylayer = TRUE, zerocol = FALSE) {
    entry <- list(...)
    ## Are any entries lists themselves? If so, make them all lists.
    listentries <- any(sapply(entry, is.list))
    if (listentries) entry <- lapply(entry, as.list)
    ## Get maximum ncol or nrow:
    lens <- sapply(entry, length)
    maxind <- max(lens)
    ## What should go in the empty part of the matrix?
    zero <- if (listentries) list(NULL) else 0
    #### Are any entries strings? If so, empty part is "".
    charentries <- any(sapply(entry, is.character))
    if (charentries) zero <- ""
    ## Append "zeroes" and construct matrix.
    numzeroes <- maxind - lens
    if (bylayer) {
        for (i in 1:length(entry)) {
            entry[[i]] <- c(entry[[i]], rep(zero, numzeroes[i]))
        }

    } else {
        for (i in 1:length(entry)) {
            entry[[i]] <- c(rep(zero, numzeroes[i]), entry[[i]])
        }
    }
    res <- sapply(entry, identity)
    if (!bylayer) res <- t(res)
    if (zerocol) res <- cbind(array(zero, c(nrow(res), 1)), res)
    res
}
