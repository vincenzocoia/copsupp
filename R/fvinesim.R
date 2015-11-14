#' Friendly vine simulation
#'
#' A (hopefully) user-friendly function to simulate from a vine copula.
#' Essentially a wrapper
#' for \code{\link{rvinesimvec2}}. (**Doesn't yet work for truncated vines,
#' because I can't figure out how to do truncated vines using
#' \code{\link{rvinesimvec2}})
#'
#' @param n Number of observations to generate
#' @param A Vine array matrix, possibly truncated.
#' @param cops
#' Vector of strings of the copula names for each edge, in "reading
#' order" corresponding to \code{A} from top-left to bottom-right.
#' To have the same copula for the entire tree, just name the copula (i.e.
#' length 1). To have the same copula for each tree, just specify the copulas
#' for each tree in that order (of length = number of trees).
#' @param cpars List, where each entry is a vector of parameters
#' corresponding to the copula in \code{cops}. Optionally, if each copula
#' family only has one parameter, this could be a vector. Should have
#' \code{length(cpars)} = # of edges in the (truncated) vine.
#' @param iprint Logical, as in \code{\link{rvinesimvec2}}, which says
#' "print flag for intermediate results".
#'
#' @details
#' The vine array in \code{A} is made up of \code{diag(1:d)} on the diagonal,
#' with 0 lower triangle, and the vine array in the upper triangle; then the
#' bottom rows are optionally removed to truncate the vine. Here, \code{d} is
#' the dimension of the vine copula (\code{=ncol(A)}).
#'
#' To name the copulas, use the names as in the CopulaModel package. For
#' example, "Frank" is \code{"frk"}, and "Gumbel" is \code{"gum"}.
#' @note
#' Like in \code{\link{rvinesimvec2}}, the copula families are assumed to be
#' permutation symmetric.
#' @examples
#' ## Vine array:
#' A <- CopulaModel::Dvinearray(5)
#' ## Simulate 10 observations with Frank copulas:
#' fvinesim(10, A, cops="frk", cpars=rep(2, choose(5,2)))
#' ## Same thing, but truncate after tree 2 (doesn't work yet):
#' fvinesim(10, A[1:2, ], cops="frk", cpars=rep(2, 7))
#' @import CopulaModel
#' @export
fvinesim <- function(n, A, cops, cpars, iprint=FALSE){
  ## Get extra parameters for rvinesimvec2 function
  ## Dimension of the vine copula:
  d <- ncol(A)
  ## ntrunc: Truncation of the vine (can't be more than d-1)
  ntrunc <- min(d-1, nrow(A))
  if (ntrunc < d-1) {
    warning("In fvinesim: I haven't yet figured out how to make a truncated vine work.")
  }
  ## parvec: Vector of parameters
  parvec <- c(cpars, recursive = T)
  ## np: Dimension of the copula parameters
  np <- makeuppertri(sapply(cpars, length), ntrunc, d)
  ## qcondmat and pcondmat:
  #### How many copulas should there be?
  numcops <- choose(d, 2) - choose(d - ntrunc, 2)
  #### Construct vector of copulas of that length.
  if (length(cops) == 1) {  # One copula given. Applies to whole vine.
    cops <- rep(cops, numcops)
  }
  if (length(cops) == ntrunc) {
    comp <- character(0)
    for (tree in 1:ntrunc) {    # One copula per tree given.
      comp <- c(comp, rep(cops[tree], d-tree))
    }
    cops <- comp
  }
  #### Now make the desired matrices:
  qmat <- makeuppertri(paste0("qcond", cops), ntrunc, d, blanks="")
  pmat <- makeuppertri(paste0("pcond", cops), ntrunc, d, blanks="")
  ## Input arguments into CopulaModel function:
  rvinesimvec2(n, A, ntrunc=ntrunc, parvec=parvec, np=np,
               qcondmat=qmat, pcondmat=pmat, iprint=iprint)
}
