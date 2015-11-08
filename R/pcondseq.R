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
#' },
#'
#' It should return a vector.
#' @export
pcondseq <- function(ord, xdat, fX, Fcond = NULL){
    if (is.vector(xdat)) xdat <- matrix(xdat, nrow = 1)
    p <- length(ord)
    ## Permute xdat so that the variables are in order of 'ord'.
    xdat <- xdat[, ord]
    ## Change fX so that it accepts a permuted vector.
    ## Note: Will need to re-permute the permuted vector back to normal. So get
    ##        the inverse permutation first.
    ordinv <- sapply(1:p, function(i) which(ord == i))
    fXperm <- function(xvecperm) {
        xvec <- xvecperm[ordinv]
        fX(xvec)
    }
    ## Find the marginal distributions of 1, 1:2, ..., 1:p by integrating.
    ##  We might not need all of them, but the integration only happens
    ##  when the function is called anyway.
    ##  (Yes we're building the list backwards, because I can't fill in the
    ##  list starting at the end. Just reverse it after.)
    fXseq <- list(fXperm)
    for (i in 1 + seq_len(p-1)) {  # is 2:p as long as p>1; integer(0) otherwise
        fXseq[[i]] <- function(xvec) {
            g <- function(xlast) fXseq[[i-1]](c(xvec, xlast))
            g <- Vectorize(g)
            integrate(g, -Inf, Inf)$value
        }
    }
    fXseq <- rev(fXseq)
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

# pcondX.vine <- function(ord, xdat, rvinefit) {
#
# }
