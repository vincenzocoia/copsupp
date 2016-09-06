#' Compute Conditional Probabilities Sequentially
#'
#' Given a vine describing the joint distribution of some uniform data, this
#' function sequentially evaluates conditional distributions given previous
#' variables in a specified order. See details.
#'
#' @param dat Data matrix with uniform scores.
#' @param ord order of variables
#' @param rv A regular vine, supposedly describing the joint
#' distribution of \code{dat}.
#' @details Suppose \code{Ui} denotes random variable \code{i}.
#' If you specify \code{ord=c(5, 2, 3, 1)}, you'll get the cdfs of
#' \code{U5}; \code{U2|U5}; \code{U3|(U2,U5)}; and \code{U1|(U3,U2,U5)},
#' evaluated at the data.
#'
#' In general, this function computes:
#'
#' \itemize{
#'      \item \code{U[ord[1]]}
#'      \item \code{U[ord[2]]|U[ord[1]]}
#'      \item \code{U[ord[3]]|U[ord[1:2]]}
#'      \item ...
#' }
#'
#' @return Matrix with rows being observations in \code{dat}, and columns being
#' \code{U[ord[1]]}; \code{U[ord[2]]|U[ord[1]]};
#' \code{U[ord[3]]|U[ord[1:2]]}; etc.
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
    if (length(ord) == 0) return(matrix(ncol=0, nrow=nrow(dat)))
    if (!all(ord %in% rv$G[1, ]))
        stop("Some variables are not described by the 'rvine' object.")
    sapply(1:length(ord), function(i) {
        pcondrvine(dat, rv, var = ord[i], condset = ord[seq_len(i-1)])
    })
}
