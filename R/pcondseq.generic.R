#' Find sequential conditional cdfs -- From density
#'
#' From a joint density, evaluates the sequential conditional cdfs
#' of a selection of variables. For example, if you choose variables 4, 2, 6, 5
#' (in that order), then it evaluates the cdfs of 4, 2|4, 6|(4,2), and 5|(4,2,6).
#'
#' @param ord Integer vector; variables in the order of finding conditional cdfs
#' (so we'll find \code{ord[1]}, \code{ord[2]|ord[1]}, etc.).
#' @param dat \code{p}-length vector, or \code{p}-columns matrix of predictors.
#' @param fX Function; the density of the covariates. Should accept a vector
#' with each component in the space (-Inf, Inf) and return a non-negative real.
#' @param Fcond If you already have some of the conditional distributions,
#' put them here in a list to speed up the algorithm. Include them either as a
#' function (see details), or a vector already
#' evaluated at the data. Make a \code{NULL} entry if you don't have that cdf.
#' The \code{k}th entry should correspond to the cdf of
#' \code{ord[k]|ord[1:(k-1)]}.
#' @return If \code{dat} is a vector, returns a
#' vector of evaluated cdfs of predictors
#' \code{ord[1]}, \code{ord[2]|ord[1]}, ..., \code{ord[p]|ord[1:(p-1)]}.
#'
#' If \code{dat} is a matrix, returns a matrix of such evaluated cdfs.
#' @note If some of your covariates don't have support on (-Inf, Inf), be sure
#' that the density still evaluates properly (to zero) outside of the support,
#' because this function integrates from -Inf to Inf.
#'
#' This function is intended as a preliminary step before connecting a response
#' to predictors in some order.
#' @details If you include a function as an entry in \code{Fcond}, it should
#' accept a vector representing variables \code{ord[c(k, 1:(k-1))]}.
#' @examples
#' (dat <- matrix(rnorm(3*5), ncol = 3))
#' pdf <- function(x) prod(dnorm(x))
#' pcondseq.generic(c(3, 1), dat, fX=pdf, Fcond = list(pnorm, NULL))
#' @seealso \code{\link{pcondseq.vine}}
#' @export
pcondseq.generic <- function(ord, dat, fX, Fcond = NULL){
    if (is.vector(dat)) dat <- matrix(dat, nrow = 1)
    p <- length(ord)  # Could be less than ncol(dat)
    pdat <- ncol(dat)
    ## Permute dat so that the variables are in order of 'ord'.
    dat <- dat[, ord]
    if (p == 1) dat <- matrix(dat, ncol = 1)
    if (is.vector(dat)) dat <- matrix(dat, nrow = 1)
    ## Change fX so that it accepts a vector in the order of ord.
    notused <- setdiff(1:pdat, ord)
    if (length(notused) > 0) {
        fXperm <- function(xperm) {
            integrand <- function(xnotused) {
                xfull <- rep(NA, pdat)
                xfull[ord] <- xperm
                xfull[notused] <- xnotused
                fX(xfull)
            }
            integrate.mv(integrand, rep(-Inf, pdat-p), rep(Inf, pdat-p))
        }
    } else {
        ## Inverse order:
        ordinv <- sapply(1:p, function(i) which(ord == i))
        fXperm <- function(xperm) fX(xperm[ordinv])
    }
    ## Work with one observation at a time.
    res <- list()
    for (i in 1:nrow(dat)) {
        res[[i]] <- numeric(0)
        xvec <- dat[i, ]
        ## Get integrands needed to compute F1, F2|1, F3|1:2, ..., Fp|1:(p-1)
        integrand <- list()
        for (j in 1:p) {
            ## Get the integrand needed to compute Fj|1:(j-1):
            ## integrate-out the upper variables, and evaluate at the lower ones
            integrand[[j]] <- function(xj) {
                pareval <- function(xup = numeric(0))
                    fXperm(c(xvec[seq_len(j-1)], xj, xup))
                integrate.mv(pareval, rep(-Inf, p-j), rep(Inf, p-j))
            }
            ## Integrate to find the Fcond, unless it's already provided.
            if (is.null(Fcond[[j]])) {
                res[[i]][j] <- integrate.mv(integrand[[j]], -Inf, xvec[j]) /
                    integrate.mv(integrand[[j]], -Inf, Inf)
            } else {
                if (is.function(Fcond[[j]])) {
                    res[[i]][j] <- Fcond[[j]](c(xvec[j], xvec[seq_len(j-1)]))
                } else {
                    res[[i]][j] <- Fcond[[j]][i]
                }
            }
        }
    }
    do.call(rbind, res)
}
