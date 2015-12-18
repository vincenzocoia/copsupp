#' Simulate from a Regular Vine
#'
#' Simulate from a Regular Vine. Essentially a wrapper
#' for \code{\link{rvinesimvec2}}. NOTE: Marginals not working yet, so
#' data are marginally Uniform.
#'
#' @param n Number of observations to generate
#' @param rv Regular vine object, complete.
#' @param iprint Logical, as in \code{\link{rvinesimvec2}}, which says
#' "print flag for intermediate results".
#' @details
#' The vine variables \code{vars(rv)} are the column numbers of the outputted
#' data. So \code{NA} columns are introduced if necessary. That is,
#' if the variables are outside of the set \code{{1:length(vars(rv))}},
#' then the columns of the outputted data outside of the set
#' \code{{1:max(vars(rv))}} will be \code{NA}.
#' @note
#' Like in \code{\link{rvinesimvec2}}, the copula families are assumed to be
#' permutation symmetric.
#' @return
#' Matrix with \code{n} rows and \code{max(vars(rv))} columns.
#' @examples
#' rv <- rvine(CopulaModel::Dvinearray(5), "frk", makeuppertri(2, 4, 5))
#' set.seed(123)
#'
#' ## Simulate 10 observations:
#' rrvine(10, rv)
#'
#' ## Simulate 0 and 1 variables too:
#' rrvine(10, subset(rv, integer(0)))
#' rrvine(10, subset(rv, 1))
#'
#' ## Vine variables outside of the set 1,2,3,...?
#' rrvine(10, relabel(rv, c(3, 8, 1, 6, 2)))
#' @import CopulaModel
#' @export
rrvine <- function(n, rv, iprint=FALSE){
    A <- rv$A
    copmat <- rv$copmat
    cparmat <- rv$cparmat
    v <- vars(rv)
  ## Get extra parameters for rvinesimvec2 function
  ## Dimension of the vine copula:
  d <- ncol(A)
  if (d == 0) return(matrix(ncol = 0, nrow = n))
  if (d == 1) {
      res <- matrix(NA, nrow = n, ncol = v)
      res[, v] <- runif(n)
      return(res)
  }
  if (d == 2) {
      res <- matrix(NA, nrow = n, ncol = max(v))
      u1 <- runif(n)
      res[, v[1]] <- u1
      qcond <- get(paste0("qcond", copmat[1, 2]))
      u2 <- qcond(runif(n), u1, cparmat[1, 2][[1]])
      res[, v[2]] <- u2
      return(res)
  }
  ## ntrunc: Truncation of the vine (can't be more than d-1)
  ntrunc <- nrow(A) - 1
  ## If there are any independence copulas, trick rvinesimvec2() by putting
  ##  a copula family with parameter that gives independence copula.
  for (i in 1:ntrunc) for (j in (i+1):d) {
      if (copmat[i, j] == "indepcop") {
          copmat[i, j] <- "new"
          cparmat[i, j][[1]] <- 0
      }
  }
  ## np: Dimension of the copula parameters
  np <- apply(cparmat, 1:2, length)
  ## qcondmat, pcondmat, parvec:
  #### How many copulas should there be?
  numcops <- choose(d, 2) - choose(d - ntrunc, 2)
  #### parvec: Vector of parameters
  parvec <- c(t(cparmat)[lower.tri(t(cparmat))], recursive = TRUE)
  #### Now make the desired matrices:
  qmat <- apply(copmat, 1:2, function(cop) paste0("qcond", cop))
  pmat <- apply(copmat, 1:2, function(cop) paste0("pcond", cop))
  qmat[!upper.tri(qmat)] <- ""
  pmat[!upper.tri(pmat)] <- ""
  ## Relabel A so that the variable name is the variable order.
  A <- relabel(rv)$A
  ## Inflate A so that it's dxd:
  if (ntrunc < d-1) {
      A <- rbind(A, matrix(0, nrow = d-(ntrunc+1), ncol = d))
      diag(A) <- 1:d # Don't worry about removing A[ntrunc+1, (ntrunc+1):d].
  }
  ## Input arguments into CopulaModel function:
  res <- rvinesimvec2(n, A, ntrunc=ntrunc, parvec=parvec, np=np,
               qcondmat=qmat, pcondmat=pmat, iprint=iprint)
  ## The user wanted the variables to be printed in the order of v:
  res2 <- matrix(NA, nrow = n, ncol = max(v))
  res2[, v] <- res
  res2
}
