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
#' @param xdat \code{p}-length vector, or \code{p}-columns matrix of predictors.
#' @param fX Function; the density of the covariates. Should accept a vector
#' with each component in the space (-Inf, Inf) and return a non-negative real.
#' @param Fcond If you already have some of the conditional distributions,
#' put them here in a list to speed up algorithm (either as a
#' function (see details), or a vector already
#' evaluated at the data). Make a \code{NULL} entry if you don't have that cdf.
#' The \code{k}th entry should correspond to the cdf of
#' \code{ord[k]|ord[1:(k-1)]}.
#' @return If \code{xdat} is a vector, returns a
#' vector of evaluated cdfs of predictors
#' \code{ord[1]}, \code{ord[2]|ord[1]}, ..., \code{ord[p]|ord[1:(p-1)]}.
#'
#' If \code{xdat} is a matrix, returns a matrix of such evaluated cdfs.
#' @note If some of your covariates don't have support on (-Inf, Inf), be sure
#' that the density still evaluates properly (to zero) outside of the support,
#' because this function integrates from -Inf to Inf.
#' @details If you include a function as an entry in \code{Fcond}, it should
#' either:
#'
#'  \enumerate{
#'      \item accept a vector if it's the cdf of a
#'      single variable (that is \code{ord[1]}), or
#'      \item accept a matrix if it's the cdf of \code{ord[k]|ord[1:(k-1)]}, with
#'      columns representing variables \code{ord[c(k, 1:(k-1))]}.
#' }
#'
#' It should return a vector.
#' @export
pcondseq.generic <- function(ord, xdat, fX, Fcond = NULL){
    if (is.vector(xdat)) xdat <- matrix(xdat, nrow = 1)
    p <- length(ord)  # Could be less than ncol(xdat)
    pdat <- ncol(xdat)
    ## Permute xdat so that the variables are in order of 'ord'.
    xdat <- xdat[, ord]
    if (p == 1) xdat <- matrix(xdat, ncol = 1)
    if (is.vector(xdat)) xdat <- matrix(xdat, nrow = 1)
    ## Change fX so that it accepts a permuted vector.
    ## Note: Will need to re-permute the permuted vector back to normal. So get
    ##        the inverse permutation first.
    ordinv <- sapply(1:p, function(i) which(ord == i))
    if (p < pdat) {
        ## Need to integrate-out unused variables.
        integration <- list()
        fXperm <- function(xvecperm) {
            for (int in setdiff(1:pdat, ord)) {
                integration <- c(integration, )
            }
            fX(xvecperm[ordinv])
        }
    } else {
        fXperm <- function(xvecperm) fX(xvecperm[ordinv])
    }
    ## Work with one observation at a time.
    for (i in 1:nrow(xdat)) {
        xvec <- xdat[i, ]
        ## Get integrands needed to compute Fp|1:(p-1), F(p-1)|1:(p-2), etc.
        integrand <- list()
        pareval <- list() # partially evaluated pdf
        for (j in 1:(p-1)) {
            ## Evaluate at conditional variables:
            pareval[[j]] <- function(xj, xupper) fXperm(c(xvec[seq_len(j-1)], xj, xupper))
            ## Integrate-out the upper variables, and get the (vectorized)
            ##  integrand in variable j.
            integrand <- function(x1j) fXperm(c(x1j, xvec[seq_len()]))
            for (int in j+seq_len(p-j)) { # Variables (j+1):p need integrating

            }
        }
        pareval[[p]] <- function(xp) fXperm(c(xvec[seq_len(p-1)], xp))

    }


    ## Find the marginal distributions of 1, 1:2, ..., 1:p by integrating.
    ##  We might not need all of them, but the integration only happens
    ##  when the function is called anyway.
    ##  (Yes we're building the list backwards, because I can't fill in the
    ##  list starting at the end. Just reverse it after.)
    fXseqrev <- list(fXperm)
    for (i in 1 + seq_len(p-1)) {  # is 2:p as long as p>1; integer(0) otherwise
        fXseqrev[[i]] <- function(xvec) {
            g <- function(xlast) fXseqrev[[i-1]](c(xvec, xlast))
            gg <- Vectorize(g)
            integrate(gg, -Inf, Inf)$value
        }
    }
    fXseq <- rev(fXseqrev)
    ## Find the conditional cdfs one-by-one:
    sapply(1:p, function(k){
        ## We're on the cdf of k|1:(k-1). Get conditioned variable indices:
        cond <- ord[seq_len(k-1)] # Could be empty.
        ## Get evaluated conditional distribution:
        ## If this conditional cdf is already given, no need for integration.
        if (!is.null(Fcond[[k]])) {
            ## Are the cdfs already evaluated?
            if (is.vector(Fcond[[k]])) {
                ## Yes. Just use them.
                res <- Fcond[[k]]
            } else {
                ## No. A vectorized function was entered, so just evaluate.
                res <- Fcond[[k]](xdat[, c(k, cond)])
            }
        } else {
            ## Need to integrate. So, need to work with one observation at a time.
            res <- sapply(1:nrow(xdat), function(row){
                ## Make density of variables c(cond, k), as a function of
                ##  variable k (i.e. evaluated at cond variables)
                xvarcond <- xdat[row, cond]
                xvark <- xdat[row, k]
                fk <- function(xk_) fXseq[[k]](c(xvarcond, xk_))
                ## Integrate to evaluate conditional cdf.
                fk <- Vectorize(fk)
                integ <- integrate(fk, -Inf, xvark)
                area <- integrate(fk, -Inf, Inf)
                if (integ$message != "OK")
                    stop (paste0("Integrate error, variable ", k, ": ", integ$message))
                integ$value / area$value
            })
        }
        res
    })
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
