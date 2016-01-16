#' Quantile Function for Y|X
#'
#' Computes the conditional quantile function of a response
#' at specified quantile indices. The model is obtained by specifying
#' copulas that link together the response with a sequence of
#' "partial" predictors (conditional on the previous predictors). Intended
#' for internal use.
#'
#' @param tau Vector of quantile indices to evaluate at. Or,
#' different vectors (of the same length) can be specified
#' for each observation by making this a matrix.
#' @param ucond Matrix of uniform "partial" predictors. Columns are ordered
#' by linkage with the response; rows are observations.
#' @param cops Vector of copula model names, corresponding to the columns of
#' \code{ucond}.
#' @param cpars List of parameter vectors corresponding to \code{cops}.
#' @param QYfitted Vectorized fitted (marginal) quantile
#' function of the response.
#' @return A matrix, with columns being the quantile levels, and rows
#' corresponding to the observations in \code{ucond}.
#' @import CopulaModel
#' @examples
#' (dat <- matrix(runif(3*6), ncol = 3))
#' tau <- c(0.1, 0.3, 0.7, 0.9)
#' cops <- c("gum", "bvtcop", "frk")
#' cpars <- list(3.5, c(0.6, 3), 2.5)
#' QYgX(tau, dat, cops = cops, cpars = cpars, QYfitted = qexp)
#' @export
QYgX <- function(tau, ucond, cops, cpars, QYfitted) {
    ## We'll evaluate the quantile function recursively. Starting
    ##  at the "tail" of ucond, we'll gradually "modify" tau
    ##  according to the copula models.
    n <- nrow(ucond)
    d <- ncol(ucond)
    if (is.vector(tau)) tau <- matrix(tau, nrow = n, ncol = length(tau), byrow = TRUE)
    if (d == 0) {
        ## There are no more predictors, so there's no need to modify tau.
        ## Evaluate at the marginal.
        return(apply(tau, 2, QYfitted))
    }
    ## There are still some predictors. Use the last one to modify tau.
    u <- ucond[, d]
    qcond <- get(paste0("qcond", cops[d]))
    qcondfit <- function(p) qcond(p, u, cpars[[d]])
    taunew <- apply(tau, 2, qcondfit)
    ## Remove a data column and evaluate the quantile function at taunew
    ucondnew <- ucond[, -d]
    if (d == 2) ucondnew <- matrix(ucondnew, ncol = 1)
    return(QYgX(taunew, ucondnew, cops, cpars, QYfitted))
}
