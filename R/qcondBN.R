#' Conditional Distribution: Bayesian Network Method
#'
#' Finds the quantiles of the distribution of a response given predictors,
#' when the response is linked to the predictors via a
#' Bayesian Network (see details), and the joint distribution of the
#' predictors is specified.
#'
#' @param tau Vector of quantile indices to evaluate at.
#' @param cop Vector of copula family names that link the response with the
#' predictors. Entry \code{i} links \code{(Y, X[i])|X[1:(i-1)]}.
#' @param cpar List of copula parameters corresponding to \code{cop}. Or,
#' vector of parameters if each copula family in \code{cop} has one parameter.
#' @param Fcond Matrix of the evaluated conditional cdfs of the predictors,
#' or a vector if there's only one observation.
#' Rows correspond to observations, and column \code{i} should contain evaluated
#' cdfs of \code{X[i]|X[1:(i-1)]}.
#' @param QY Marginal quantile function for the response.
#' @return Matrix of conditional quantiles of the response given the
#' predictors, or a vector if there's only one observation in \code{Fcond}.
#' Columns correspond to quantile indices \code{tau}.
#' @details
#' The conditional distribution of response \code{Y} conditional on \code{p}
#' predictors \code{X} is built by specifying the pairwise distributions of
#' \code{Y} with the responses:
#'
#' \enumerate{
#'      \item \code{(Y,X1)}
#'      \item \code{(Y,X2)|(X1)}
#'      \item \code{(Y,X3)|(X1,X2)}
#'      \item \code{(Y,X4)|(X1,X2,X3)}
#' }
#'
#' For each, a copula should describe the dependence.
#'
#' @examples
#' ## Setup predictor distribution: D-vine
#' set.seed(123)
#' A <- truncvarray(CopulaModel::Dvinearray(4), 2)
#' copmat <- makeuppertri("frk", 2, 4, "")
#' cparmat <- makeuppertri(5:1, 2, 4)
#' dat <- fvinesim(10, A, copmat, cparmat)
#'
#' ## Link Y with X1,...,X4 in this order:
#' ord <- c(3, 4, 2, 1)
#' Fcond <- pcondseq.vine(ord, dat,
#'                        rvinefit=list(A=A, copmat=copmat, cparmat=cparmat))
#'
#' ## ...and with these copulas:
#' Ycop <- c("frk", "gum", "mtcj", "bvncop")
#' Ypar <- c(4, 2, 2, 0.7)
#'
#' ## Find 0.9-, 0.95-, and 0.99-quantiles of Y|X at the data, with Y~Exp(1):
#' qcondBN(c(0.9, 0.95, 0.99), cop=Ycop, cpar=Ypar, Fcond=Fcond, QY=qexp)
#'
#' ## Try one quantile index:
#' qcondBN(0.9, cop=Ycop, cpar=Ypar, Fcond=Fcond, QY=qexp)
#'
#' ## Try one observation:
#' qcondBN(c(0.9, 0.95, 0.99), cop=Ycop, cpar=Ypar, Fcond=Fcond[1, ], QY=qexp)
#' @seealso See \code{\link{pcondseq.vine}} and \code{\link{pcondseq.generic}}
#' for computing \code{Fcond}.
#' @export
qcondBN <- function(tau, cop, cpar, Fcond, QY = identity) {
    if (!is.matrix(Fcond)) return(qcondBN.oneobs(tau, cop, cpar, Fcond, QY))
    p <- ncol(Fcond)
    n <- nrow(Fcond)
    K <- length(tau)
    if (p == 0) return(matrix(QY(tau), ncol = K, nrow = n, byrow = TRUE))
    res <- apply(Fcond, 1, function(row) qcondBN.oneobs(tau, cop, cpar, row, QY))
    if (K == 1) return(matrix(res, ncol = 1)) else return(t(res))
}

qcondBN.oneobs <- function(tau, cop, cpar, Fcondvec, QY = identity) {
    p <- length(Fcondvec)
    K <- length(tau)
    if (p == 0) return(QY(tau))
    ## Recursive form of conditional qf of response:
    thiscop <- get(paste0("qcond", cop[p]))
    if (is.list(cpar)) {
        thiscpar <- cpar[[p]]
    } else {
        thiscpar <- cpar[p]
    }
    newtau <- thiscop(tau, Fcondvec[p], thiscpar)
    qcondBN.oneobs(newtau, cop, cpar, Fcondvec[-p], QY = QY)
}
