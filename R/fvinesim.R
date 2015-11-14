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
#' @param cops Upper-triangular matrix of copula names (like "frk" or "gum"),
#' corresponding to the edges in vine array \code{A}. Or, a single name
#' if it applies to all edges in \code{A}, or a vector of length
#' \code{nrow(A)-1} where entry \code{i} corresponds to edges in row \code{i}
#' of \code{A}.
#' @param cpars Matrix of copula parameters, with entries corresponding to
#' entries in \code{cops}. Optionally, if only one copula is being fit for
#' the entire vine, could be the copula parameter (a vector).
#' @param iprint Logical, as in \code{\link{rvinesimvec2}}, which says
#' "print flag for intermediate results".
#'
#' @details
#' To truncate a vine array, use \code{\link{trunc.varray}}.
#'
#' To make a matrix of copulas, use \code{\link{makeuppertri}}. Some copula
#' parameters may have more or less parameters than 1 -- in this case,
#' enlist the help of \code{\link{makeuppertri.list}}.
#'
#' To name the copulas, use the names as in the CopulaModel package. For
#' example, "Frank" is \code{"frk"}, and "Gumbel" is \code{"gum"}.
#' @note
#' Like in \code{\link{rvinesimvec2}}, the copula families are assumed to be
#' permutation symmetric.
#' @examples
#' ## Vine array:
#' A <- CopulaModel::Dvinearray(5)
#' A <- relabel.varray(A, c(3, 5, 1, 2, 4))
#' ## Simulate 10 observations with Frank copulas:
#' set.seed(123)
#' fvinesim(10, A, cops="frk", cpars=2)
#' ## Same thing, but 2-truncated:
#' A <- trunc.varray(A, 2)
#' set.seed(123)
#' fvinesim(10, A, cops="frk", cpars=2)
#' ## Notice that variables 3,5,1 -- the first three generated --- are the same
#' ##  as the complete vine, since they are only linked by 2 trees anyways.
#' @import CopulaModel
#' @export
fvinesim <- function(n, A, cops, cpars, iprint=FALSE){
  ## Get extra parameters for rvinesimvec2 function
  ## Dimension of the vine copula:
  d <- ncol(A)
  ## ntrunc: Truncation of the vine (can't be more than d-1)
  ntrunc <- nrow(A) - 1
  ## np: Dimension of the copula parameters
  if (is.vector(cpars)) {
      np <- makeuppertri(length(cpars), ntrunc, d)
  } else {
      if (is.list(cpars[1,1])){
          np <- apply(cpars, 1:2, function(t) length(t[[1]]))
      } else {
          np <- apply(cpars, 1:2, length)
      }
  }
  ## qcondmat, pcondmat, parvec:
  #### How many copulas should there be?
  numcops <- choose(d, 2) - choose(d - ntrunc, 2)
  #### parvec: Vector of parameters
  parvec <- c(t(cpars), recursive = T)
  if (is.vector(cpars)) parvec <- rep(parvec, numcops)
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
  ## Relabel A so that the variable name is the variable order.
  labs_orig <- varray.vars(A)
  A <- relabel.varray(A)
  ## Inflate A so that it's dxd:
  if (ntrunc < d-1) {
      A <- rbind(A, matrix(0, nrow = d-ntrunc, ncol = d))
      diag(A) <- 1:d # Don't worry about removing A[ntrunc+1, (ntrunc+1):d].
  }
  ## Input arguments into CopulaModel function:
  res <- rvinesimvec2(n, A, ntrunc=ntrunc, parvec=parvec, np=np,
               qcondmat=qmat, pcondmat=pmat, iprint=iprint)
  ## The user wanted the variables to be printed in the order of labs_orig:
  res[, invert.perm(labs_orig)]
}
