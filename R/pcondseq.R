#' Find sequential conditional cdfs
#'
#' From a \code{p}-variate joint distribution,
#' finds the conditional cdfs of variables \code{1}, \code{2|1}, ...,
#' \code{p|1:(p-1)} evaluated at some data.
#' This function is intended as a preliminary step before connecting a response
#' to covariate \code{1}, then \code{2|1}, ...,\code{p|1:(p-1)}.
#' Allows for permutations of \code{1:p} too.
#'
#' @param ord Integer vector; variables in the order that you'll be linking them
#' up with the response (so we'll find \code{ord[1]}, \code{ord[2]|ord[1]}, etc.).
#' @param dat \code{p}-length vector, or \code{p}-columns matrix of predictors.
#' @param fX Function; the density of the covariates. Should accept a vector
#' with each component in the space (-Inf, Inf) and return a non-negative real.
#' @param Fcond If you already have some of the conditional distributions,
#' put them here in a list to speed up algorithm (either as a
#' function (see details), or a vector already
#' evaluated at the data). Make a \code{NULL} entry if you don't have that cdf.
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
#' @details If you include a function as an entry in \code{Fcond}, it should
#' accept a vector representing variables \code{ord[c(k, 1:(k-1))]}.
#' @examples
#' dat <- c(0.4, 0.5, 0.6)
#' pdf <- function(x) as.integer(all(x>0) & all(x<1))
#' pcondseq.generic(c(3, 1), dat, fX=pdf, Fcond = list(identity, NULL, NULL))
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
    res <- list(numeric(0))
    for (i in 1:nrow(dat)) {
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

#' Temp
#'
#' @param ord Integer vector; variables in the order that you'll be linking them
#' up with the response (so we'll find \code{ord[1]}, \code{ord[2]|ord[1]}, etc.).
#' @param xdat Vector of a single observation, or matrix of multiple observations
#' of the variables.
#' @param rvinefit The vine fit to data \code{xdat}. See details.
#' @param FX List of vectorized functions that are the marginals corresponding
#' to the columns of \code{xdat}. Or just one function if it's a common function.
#' Default is \code{identity} so that \code{xdat} can be uniform data if you want.
#' @details The argument \code{rvinefit} can either be:
#'
#' \enumerate{
#'      \item the output of \code{\link{VineCopula::RVineCopSelect}} (Version 1.6), or
#'      \item a list with the following named entries:
#'      \itemize{
#'          \item \code{A}: The vine array, as used in the
#'          package \code{\link{CopulaModel}}.
#'          \item \code{copmat}: An upper-triangular matrix of names of the
#'          bivariate copula models used in the vine.
#'          \item \code{cparmat}: An upper-triangular matrix of copula parameters
#'          to use in the corresponding copula model in \code{rvinefit$copmat}.
#'          Each entry should be a vector with length = the number of parameters
#'          for that copula model. See \code{\link{makeuppertri.list}} for help.
#'      }
#' }
pcondseq.vine <- function(ord, xdat, rvinefit, FX = identity) {
    ## Standardize input
    if (is.vector(xdat)) xdat <- matrix(xdat, ncol = length(xdat))
    if (nrow(xdat) == 0) if (ncol(xdat) == 0) return(numeric(0)) else return(xdat)
    p <- length(ord)  # May be < ncol(xdat).
    if (length(FX) == 1) FX <- rep(list(FX), p)
    ## re-order xdat and marginals so that they're in the order of ord.
    xdat <- xdat[, ord]
    if (is.vector(xdat)) xdat <- matrix(xdat, ncol = 1)
    FX <- FX[ord]
    ## Uniformize data
    udat <- xdat
    for (col in 1:p) udat[, col] <- FX[[col]](xdat[, col])
    ## Don't go any further if p==1
    if (p == 1) return(udat)
    ## Vine-specific input:
    A <- rvinefit$A
    ## Get sub-vine arrays for ord[1], ..., ord[k], if they exist.
    subA <- list()
    isleaf <- logical(0)
    isD <- logical(0)
    for (k in 1:p) {
        subA[[k]] <- rvinesubset(A, ord[1:k])
        ## Which have ord[k] as a leaf?
        isleaf <- c(isleaf, is.vineleaf(subA[[k]], ord[[k]]))
        ## Which are D-Vines? that have ord[k] as a leaf?
        isD <- c(isD, is.dvine(subA[[k]]))
    }
    ## Find conditional cdfs (I'll have to construct the list backwards first):
    Fcondlist <- list()
    for (k in p:2) {
        subAk <- subA[[k]]
        if (is.null(subAk)) {
            ## Need to integrate the density.
        } else {
            if (isleaf[k] & isD[k]) {
                ## Can use the formula for D-Vine -- no integration.

            } else {
                ## Can integrate the conditional density.
            }
        }
    }
    #### And finally, for k == 1,
    Fcondlist <- list(udat[, 1])
    ## Reverse, and collapse into a matrix.
    Fcondlist <- rev(Fcondlist)
    do.call(cbind, Fcondlist)
}
