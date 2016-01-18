#' Find sequential predictors
#'
#' For internal use.
#'
#' @param dat Data matrix with uniform scores.
#' @param ord order of variables
#' @param rv A regular vine, supposedly fitted to \code{dat}.
#' @return Matrix with rows being observations, and columns being
#' X[ord[1]]; X[ord[2]]|X[ord[1]]; X[ord[3]]|X[ord[1:2]]; etc.
#' @note Expecting smart input: ord must be at least of length 1, with
#' variables in the specified rvine.
#' @examples
#' ## Make a vine
#' G <- AtoG(CopulaModel::Dvinearray(4))
#' rv <- rvine(G, "frk", 3)
#'
#' ## Generate sample
#' set.seed(123)
#' dat <- rrvine(10, rv)
#'
#' ## Try it out:
#' pcondseq(dat, c(4, 3, 2), rv)
#' pcondseq(dat, 1, rv)
#' @export
pcondseq <- function(dat, ord, rv) {
    sapply(1:length(ord), function(i) {
        pcondrvine(dat, rv, var = ord[i], condset = ord[seq_len(i-1)])
    })
}
