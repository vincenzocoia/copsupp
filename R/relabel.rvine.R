#' Relabel Variables in a Vine/Vine Array
#'
#' @param rv Object of type 'rvine'
#' @param G Vine array
#' @param labs Vector of new labels. The order of the labels correspond to
#' the order of the variables in the vine array \code{G}.
#' @details Similar to \code{CopulaModel::varrayperm} but allows for
#' the posibility that \code{G} is not square, as well as labels outside
#' of the set \code{{1:ncol(G)}}.
#' @return The inputted \code{rv} or \code{G} with the array matrix
#' relabelled.
#' @examples
#' (G <- AtoG(CopulaModel::Cvinearray(5))[1:3, ])
#' rv <- rvine(G, "frk", 4)
#'
#' relabelvarray(G, c(3, 2, 1, 5, 4))
#' relabel(rv, c(3, 2, 1, 5, 4))
#'
#' ## Labels don't need to be in the set {1:ncol(G)}.
#' relabelvarray(G, c(54, 234, 1, 35, 42))
#' relabel(rv, c(54, 234, 1, 35, 42))
#' @rdname relabel
#' @export
relabel.rvine <- function(rv, labs = 1:ncol(rv$G)) {
    G <- rv$G
    d <- ncol(G)
    r <- nrow(G)
    labs_orig <- vars(rv)
    ## Map original label to order
    map2order <- function(labo) which(labs_orig == labo)
    ## Map original labels to new labels
    for (row in 1:r) {
        for (col in row:d) {
            G[row, col] <- labs[map2order(G[row, col])]
        }
    }
    rvine(G, rv$copmat, rv$cparmat)
}

#' @export
relabel <- function(...) UseMethod("relabel")

#' @rdname relabel
#' @export
relabelvarray <- function(G, labs = 1:ncol(G)) {
    d <- ncol(G)
    r <- nrow(G)
    labs_orig <- G[1, ]
    ## Map original label to order
    map2order <- function(labo) which(labs_orig == labo)
    ## Map original labels to new labels
    for (row in 1:r) {
        for (col in row:d) {
            G[row, col] <- labs[map2order(G[row, col])]
        }
    }
    G
}
