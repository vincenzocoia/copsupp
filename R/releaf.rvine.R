#' "Re-leaf" a Vine Array
#'
#' Convert a vine array so that a variable(s) of your choice appears last
#' (bottom-right corner) of the vine array -- if possible.
#'
#' @param rv Regular vine object.
#' @param leaf Integer; variable you want to make a leaf of the array.
#' @return A regular vine object with \code{leaf} at the end
#' of the vine array, if such a
#' vine exists; \code{NULL} otherwise.
#' @examples
#' A <- makeuppertri(c(4,4,4,4,6,5,
#'                     6,6,6,4,4,
#'                     1,1,1,6,
#'                     5,5,1,
#'                     2,2,
#'                     3), 6, 6, incDiag=T)
#' rv <- rvine(A)
#' releaf(rv, leaf = 2)
#' releaf(rv, leaf = 1)
#' releaf(trunc(rv, 2), leaf = 1)
#'
#' ## 1:ncol(A) not required as variables.
#' rv <- subset(rv, vars(rv)[1:4])
#' releaf(rv, leaf = 1)
#' @export
releaf.rvine <- function(rv, leaf) {
    A <- rv$A
    d <- ncol(A)
    ## Empty vine case:
    if (d == 0) return(rv)
    v <- vars(rv)
    ## Stop if the leaf is not in the vine
    if (!(leaf %in% v)) stop(paste("Variable", leaf, "not found in vine."))
    ntrunc <- nrow(A) - 1
    ## If 'leaf' is already at the end of the vine, just return the vine.
    if (A[ntrunc+1, d] == leaf) return(rv)
    ## If rv is an independence vine, that's an easy case:
    if (ntrunc == 0) {
        ileaf <- which(v == leaf)
        A[1, c(ileaf, d)] <- A[1, c(d, ileaf)]
        marg <- rv$marg
        marg[c(ileaf, d)] <- marg[c(d, ileaf)]
        return(rvine(A, marg = marg))
    }
    ## Case when d=2
    if (d == 2) {
        A <- matrix(c(A[2,2], 0, A[2,2], leaf), ncol = 2)
        return(rvine(A, rv$copmat, rv$cparmat, rv$marg[2:1]))
    }
    ## Center the vine:
    rvc <- center(rv)
    ## Deal with the complete-vine case separately, in which case rvc is
    ##   in natural order.
    if (d == ntrunc + 1) {
        A <- rvc$A
        if (A[d,d] == leaf) {
            return(rvc)
        }
        if (A[d-1, d-1] == leaf) {
            Anew <- A
            copmatnew <- rvc$copmat
            cparmatnew <- rvc$cparmat
            ## Swap columns
            Anew[, d-1] <- A[, d]
            Anew[, d] <- A[, d-1]
            ## Deal with 2x2 bottom-right portion:
            Anew[d-1, d-1] <- A[d, d]
            Anew[d, d] <- A[d-1, d-1]
            Anew[d, d-1] <- 0
            Anew[d-1, d] <- Anew[d-1, d-1]
            ## Swap top d-1 portion of last two columns of copula and param matrices:
            if (!is.null(copmatnew)) {
                one <- copmatnew[1:(d-1), d-1]
                two <- copmatnew[1:(d-1), d]
                copmatnew[1:(d-1), d-1] <- two
                copmatnew[1:(d-1), d] <- one
            }
            if (!is.null(cparmatnew)) {
                one <- cparmatnew[1:(d-1), d-1]
                two <- cparmatnew[1:(d-1), d]
                cparmatnew[1:(d-1), d-1] <- two
                cparmatnew[1:(d-1), d] <- one
            }
            ## Swap marginals
            margnew <- rvc$marg
            margnew[c(d-1, d)] <- margnew[c(d, d-1)]
            return(rvine(Anew, copmatnew, cparmatnew, margnew))
        }
        return(NULL)
    }
    ## --- END special cases. ---
    Ac <- rvc$A
    ## Get variables that won't work:
    bad <- unique(as.vector(Ac[1:ntrunc, (ntrunc+1):d]))
    if (leaf %in% bad) return(NULL)
    ## Leaf if valid. Reorder vine:
    vc <- vars(rvc)
    ileaf <- which(vc == leaf)
    Ac[, c(ileaf, d)] <- Ac[, c(d, ileaf)]
    copmatc <- rvc$copmat
    copmatc[, c(ileaf, d)] <- copmatc[, c(d, ileaf)]
    cparmatc <- rvc$cparmat
    cparmatc[, c(ileaf, d)] <- cparmatc[, c(d, ileaf)]
    margc <- rvc$marg
    margc[c(ileaf, d)] <- margc[c(d, ileaf)]
    rvine(Ac, copmatc, cparmatc, margc)
}

#' @export
releaf <- function(...) UseMethod("releaf")

#' "Re-leaf" a Vine Array
#'
#' Convert a vine array so that a variable(s) of your choice appears last
#' (bottom-right corner) of the vine array -- if possible.
#'
#' @param A Vine array matrix, possibly truncated.
#' @param leaf Vector of variables that you want to
#' get new vine array for, if such a vine array exists.
#' @return A vine array (matrix) with \code{leaf} at the end, if such a
#' vine array exists; \code{NULL} otherwise. If \code{length(leaf) > 1},
#' a list of such output is returned, corresponding to the entries
#' in \code{leaf}.
#' @examples
#' A <- makeuppertri(c(4,4,4,4,6,5,
#'                     6,6,6,4,4,
#'                     1,1,1,6,
#'                     5,5,1,
#'                     2,2,
#'                     3), 6, 6, incDiag=T)
#' releaf.varray(truncvarray(A, 2))
#' releaf.varray(truncvarray(A, 3), leaf = 3:6)
#' releaf.varray(A, leaf = 5)
#'
#' ## The function doesn't require 1:ncol(A) as variables.
#' A <- truncvarray(CopulaModel::Dvinearray(5), 2)
#' A <- apply(A, 1:2, function(i) if (i==0) "" else letters[i])
#' releaf.varray(A)
#' @export
releaf.varray <- function(A, leaf=varray.vars(A)) {
    warning("'releaf.varray()' is deprecated. Please use 'releaf.rvine()' instead.")
    if (!is.matrix(A)) return(NULL)
    if (length(leaf) == 0) return(list())
    if (ncol(A) == 0) return(NULL) # Only after checking length(leaf) > 0.
    d <- ncol(A)
    ntrunc <- nrow(A) - 1
    varsset <- varray.vars(A)
    ## Case when d=2
    if (d == 2) {
        res <- lapply(leaf, function(l) {
            if (A[2,2] == l) {
                A
            } else {
                matrix(c(A[2,2], 0, A[2,2], l), ncol = 2)
            }
        })
        if (length(res) == 1) res <- res[[1]]
        return(res)
    }
    ## It's easiest to work with arrays truncated < d-1. So make d-2 the max.
    ##  After all, it's easy to go from d-2 back to the full d-1, which
    ##  we'll do later.
    putback <- FALSE
    if (ntrunc == d-1) {
        putback <- TRUE
        ntrunc <- d-2
        A <- truncvarray(A, ntrunc)
    }
    ## Centralize the vine array:
    Ac <- center.varray(A)
    leftAc <- Ac[, 1:ntrunc]
    if (!is.matrix(leftAc)) leftAc <- matrix(leftAc, ncol = ntrunc)
    ## Identify truncated region to the right of column ntrunc. They're
    ##  candidate columns/variables.
    vars <- varray.vars(Ac)
    candvars <- vars[(ntrunc+1):d]
    candbranches <- Ac[1:ntrunc, (ntrunc+1):d]
    candcols <- Ac[1:(ntrunc+1), (ntrunc+1):d]
    if (!is.matrix(candbranches))
        candbranches <- matrix(candbranches, ncol = ntrunc)
    ncols <- ncol(candcols)
    ## Variables that won't work as leaves:
    bad <- unique(as.vector(candbranches))
    works <- !(leaf %in% bad)
    ## Re-order the vine arrays:
    res <- list()
    for (i in 1:length(leaf)) {
        if (works[i]) {
            ## Leaf is valid.
            col <- which(candvars == leaf[i])
            recol <- c((1:ncols)[-col], col)
            res[[i]] <- cbind(leftAc, candcols[, recol])
            if (putback) {
                res[[i]] <- rbind(res[[i]], matrix(0, ncol = d))
                res[[i]][d,d] <- res[[i]][d-1, d]
                res[[i]][d-1,d] <- setdiff(varsset, res[[i]][, d])
            }
        } else {
            ## Leaf is invalid.
            res <- c(res, list(NULL))
        }
    }
    names(res) <- leaf
    if (length(res) == 1) res <- res[[1]]
    res
}
