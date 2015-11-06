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


#' Friendly vine simulation
#' 
#' A (hopefully) user-friendly function to simulate from a vine copula. 
#' Essentially a wrapper
#' for \code{\link{rvinesimvec2}}. (**Doesn't yet work for truncated vines,
#' because I can't figure out how to do truncated vines using 
#' \code{\link{rvinesimvec2}})
#' 
#' @param n Number of observations to generate
#' @param truncA Vine array matrix, possibly truncated (so, #columns should
#' equal vine copula dimension). 
#' @param cops Vector of strings of the copula names for each edge, in "reading
#' order" corresponding to \code{truncA} from top-left to bottom-right.
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
#' The vine array in \code{truncA} is made up of \code{diag(1:d)} on the diagonal,
#' with 0 lower triangle, and the vine array in the upper triangle; then the
#' bottom rows are optionally removed to truncate the vine. Here, \code{d} is
#' the dimension of the vine copula (\code{=ncol(truncA)}).
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
fvinesim <- function(n, truncA, cops, cpars, iprint=FALSE){
  ## Get extra parameters for rvinesimvec2 function
  ## Dimension of the vine copula:
  d <- ncol(truncA)
  ## ntrunc: Truncation of the vine (can't be more than d-1)
  ntrunc <- min(d-1, nrow(truncA))
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
  rvinesimvec2(n, truncA, ntrunc=ntrunc, parvec=parvec, np=np,
               qcondmat=qmat, pcondmat=pmat, iprint=iprint)
}