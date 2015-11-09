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

#' Subset a vine array
#'
#' Find the vine connecting a subset of variables from a bigger vine, if
#' it exists.
#'
#' @param A vine array
#' @param select Vector of variable indices in \code{diag(A)} to subset,
#' if possible.
#' @details Just a technicality:
#' by saying a subset "doesn't have an existing vine", I mean that
#' a vine can't be formed using nodes and edges from
#' the original -- not that the
#' joint distribution of the selected variables can't be created from a vine
#' (so as to say, for example, that the simplifying assumption of vines
#' doesn't hold for this distribution).
#' @return Returns a vine array of the subsetted variables, or \code{NULL} if
#' the subset don't form a vine.
#' @examples
#' A <- CopulaModel::Dvinearray(5)
#' rvinesubset(A, c(2, 3, 4))
#' rvinesubset(A, c(4, 1))
#' @export
rvinesubset <- function(A, select) {
    p <- ncol(A)
    k <- length(select)
    if (k == p) return(A)
    if (k == 1) return(matrix(select))
    diagA <- diag(A)
    notselect <- setdiff(diag(A), select)
    ## Indices to keep (i.e. what columns of A, or which of the
    ##   ordered variables to keep?)
    ikeep <- sapply(select, function(s) which(diagA == s))
    ilose <- setdiff(1:p, ikeep)
    ## Select entries to remove from the vine array (as TRUE entries).
    #### Bottom-left zeroes need removal.
    removal <- lower.tri(diag(p))
    #### "Extraneous" area needs removal (i.e. higher-level trees
    ####   that are impossible with this selection)
    removal[k:p, k:p] <- removal[k:p, k:p] | upper.tri(diag(p-k+1))
    #### Columns where our selection is not on the diagonal need removing:
    removal[, ilose] <- TRUE
    #### Select remaining variables that need removal.
    removal[, ikeep] <- removal[, ikeep] |
        apply(A[, ikeep], 1:2, function(t) t %in% notselect)
    ## There should be (k+1) choose 2 variable indices remaining if the
    ##  subsetted vine "exists".
    remain <- A[!removal]  # A vector.
    if (length(remain) == choose(k+1, 2)) {
        subA <- makeuppertri(remain, k+1, k+1, byRow = FALSE)
        subA <- subA[, -1]
        subA <- subA[-(k+1), ]
    } else {
        subA <- NULL
    }
    subA
}

#' Temp
#'
#' @param rvinefit The vine fit to data \code{xdat}. See details.
#' @details The argument \code{rvinefit} can either be:
#'
#' \enumerate{
#'      \item the output of \code{\link{VineCopula::RVineCopSelect}}, or
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
    if (is.vector(xdat)) xdat <- matrix(xdat, nrow = 1)
    p <- ncol(xdat)
    if (length(FX) != p) FX <- rep(list(FX), p)
    ## Uniformize data
    udat <- xdat
    for (col in 1:p) udat[, col] <- FX[[col]](xdat[, col])
    ## Vine-specific input:
    A <- rvinefit$A
    ## Get sub-vine arrays for ord[1], ..., ord[k], if they exist.
    ##   -Doesn't exist? Need to use generic function.
    ##   -Exists, and ord[k] isn't a leaf on a D-Vine? There's a formula for
    ##      conditional density. Just need one integral.
    ##   -Exists, and ord[k] is a leaf on a D-Vine? Closed formula for conditional cdf.
    ordinv <- sapply(1:p, function(i) which(i == ord)) # inverse permutation
    A <- varrayperm(A, ordinv[diag(A)])  # Now, diag(A) represents the order of adding k.
    subA <- list(matrix(1))
    for (k in 1+seq_len(p-2)) { # For p>=3, this is 2:(p-1).
        ## Select zeroes to remove from vine array.
        removal <- lower.tri(diag(p))
        ## Select "sensitive" triangle to remove.
        removal[k:p, k:p] <- removal[k:p, k:p] | upper.tri(diag(p-k+1))
        ## Remove unused columns
        removal[, ordinv[(k+1):p]] <- TRUE
        ## Select remaining variables that need removal.
        removal[, ordinv[1:k]] <- removal[, ordinv[1:k]] |
            apply(A[, ordinv[1:k]], 1:2, function(t) t %in% (k+1):p)
        ## There should be (k+1) choose 2 variable indices remaining if it's a
        ##  valid vine.
        remain <- A[!removal]
        if (length(remain) == choose(k+1, 2)) {
            subA[[k]] <- makeuppertri(remain, k+1, k+1, byRow = FALSE)
            subA[[k]] <- subA[[k]][, -1]
            subA[[k]] <- subA[[k]][-(k+1), ]
        } else {
            subA[[k]] <- NULL
        }
    }
    subA[[p]] <- A
    #### Which are D-Vines anyways?
    whichD <- sapply(subA, function(subA_){
        if (is.null(subA_)) return(FALSE)
        if (ncol(subA_) == 1) return(FALSE)
        ## In tree 1, nodes should appear maximum of two times.
        max(table(c(subA_[1, -1], diag(subA_)[-1]))) <= 2
    })


}
