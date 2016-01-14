#' "Re-leaf" a Vine
#'
#' Convert a vine array (\code{releafvarray}) or a vine
#' (\code{releaf.rvine}) so that a variable of your choice appears last
#' in the vine array -- if possible. This might only work for
#' regular vines, with arrays truncated in the traditional sense.
#'
#' @param G Regular vine array.
#' @param rv Object of type "rvine".
#' @param leaf Vector of variables that you want to
#' get new vine array for, if such a vine array exists.
#' @return A vine array or vine with \code{leaf} at the end, if such a
#' vine array exists; \code{NULL} otherwise.
#' @examples
#' A <- makeuppertri(c(4,4,4,4,6,5,
#'                     6,6,6,4,4,
#'                     1,1,1,6,
#'                     5,5,1,
#'                     2,2,
#'                     3), 6, 6, incDiag=T)
#' G <- AtoG(A)
#' releafvarray(G, 2)
#' releafvarray(G, 1)
#' releafvarray(truncvarray(G, 0), 1)
#' releaf(rvine(G, "frk", 3), 2)
#' @rdname releaf
#' @export
releafvarray <- function(G, leaf) {
    d <- ncol(G)
    ## Case when there are no variables
    if (d == 0) if (length(leaf) == 0) return(G) else return(NULL)
    ntrunc <- nrow(G) - 1
    ovars <- G[1, ]
    ## Case when leaf is not a variable
    if (!(leaf %in% ovars)){
        warning("Choice of leaf is not a variable in G. Returning NULL.")
        return(NULL)
    }
    ## Case when d=1, or in general, G is an independence vine array.
    if (ntrunc == 0){  # Accounts for d == 1 case too
        return(matrix(c(setdiff(ovars, leaf), leaf), nrow=1))
    }
    ## Case when d=2
    if (d == 2) {
        if (G[1,2] == leaf) {
            return(G)
        } else {
            v1 <- G[1, 2]
            v2 <- G[1, 1]
            return(makevinemat(v1, c(v2, v1)))
        }
    }
    ## We might get lucky and have 'leaf' as the end already.
    if (G[1, d] == leaf) return(G)
    ## We now have d>2. Let's center the vine array.
    Gnew <- centervarray(G)
    ## We might get lucky and have 'leaf' as the end already.
    if (Gnew[1, d] == leaf) return(Gnew)
    ## Case when vine is complete (and by now, the leaf is not at the end)
    if (ntrunc == d-1) {
        ## Since "leaf" is not at the end, it's only a valid leaf if it's
        ##  second-last.
        if (Gnew[1, d-1] == leaf) {
            neword <- Gnew[1, 1:d]
            neword[(d-1):d] <- neword[d:(d-1)]
            return(reordervarray(Gnew, neword))
        } else {
            return(NULL)
        }
    }
    ## We now have a truncated vine array with more than 2 variables.
    ##  Plus, the requested leaf is not at the end. Check whether the
    ##  requested leaf has any variables in "upstream" columns in its own column.
    ##  If so, it can't be a leaf.
    leafcol <- which(ovars == leaf)
    upstreamvars <- c(G[, (leafcol+1):d], recursive = TRUE)
    is_in_upstream <- any(sapply(G[, leafcol], function(x) x %in% upstreamvars))
    if (is_in_upstream) {
        return(NULL)
    } else {
        ## It's valid. Shift its column to the far right of the array, and
        ##  we're done.
        return(G[, c((1:d)[-leafcol], leafcol)])
    }
}

#' @rdname releaf
#' @export
releaf.rvine <- function(rv, leaf) {
    ## Re-leaf vine array first.
    G <- rv$G
    Gnew <- releafvarray(G, leaf)
    if (is.null(Gnew)) return(NULL)
    ## Transform copmat and cparmat
    copmat <- reformcopmat(rv$copmat, Gnew, G)
    cparmat <- reformcopmat(rv$cparmat, Gnew, G)
    ## Output
    rvine(Gnew, copmat, cparmat)
}

#' @export
releaf <- function(...) UseMethod("releaf")


# releaf.rvine <- function(rv, leaf) {
#     A <- GtoA(rv$G)
#     d <- ncol(A)
#     ## Empty vine case:
#     if (d == 0) return(rv)
#     v <- vars(rv)
#     ## Stop if the leaf is not in the vine
#     if (!(leaf %in% v)) stop(paste("Variable", leaf, "not found in vine."))
#     ntrunc <- nrow(A) - 1
#     ## If 'leaf' is already at the end of the vine, just return the vine.
#     if (A[ntrunc+1, d] == leaf) return(rv)
#     ## If rv is an independence vine, that's an easy case:
#     if (ntrunc == 0) {
#         ileaf <- which(v == leaf)
#         A[1, c(ileaf, d)] <- A[1, c(d, ileaf)]
#         return(rvine(AtoG(A)))
#     }
#     ## Case when d=2
#     if (d == 2) {
#         A <- matrix(c(A[2,2], 0, A[2,2], leaf), ncol = 2)
#         return(rvine(AtoG(A), rv$copmat, rv$cparmat))
#     }
#     ## Center the vine:
#     rvc <- center(rv)
#     ## Deal with the complete-vine case separately, in which case rvc is
#     ##   in natural order.
#     if (d == ntrunc + 1) {
#         A <- GtoA(rvc$G)
#         if (A[d,d] == leaf) {
#             return(rvc)
#         }
#         if (A[d-1, d-1] == leaf) {
#             Anew <- A
#             copmatnew <- rvc$copmat
#             cparmatnew <- rvc$cparmat
#             ## Swap columns
#             Anew[, d-1] <- A[, d]
#             Anew[, d] <- A[, d-1]
#             ## Deal with 2x2 bottom-right portion:
#             Anew[d-1, d-1] <- A[d, d]
#             Anew[d, d] <- A[d-1, d-1]
#             Anew[d, d-1] <- 0
#             Anew[d-1, d] <- Anew[d-1, d-1]
#             ## Swap top ntrunc-1 portion of last two columns of copula and param matrices:
#             if (!is.null(copmatnew)) {
#                 one <- copmatnew[1:(ntrunc-1), d-1]
#                 two <- copmatnew[1:(ntrunc-1), d]
#                 copmatnew[1:(ntrunc-1), d-1] <- two
#                 copmatnew[1:(ntrunc-1), d] <- one
#             }
#             if (!is.null(cparmatnew)) {
#                 one <- cparmatnew[1:(ntrunc-1), d-1]
#                 two <- cparmatnew[1:(ntrunc-1), d]
#                 cparmatnew[1:(ntrunc-1), d-1] <- two
#                 cparmatnew[1:(ntrunc-1), d] <- one
#             }
#             return(rvine(AtoG(Anew), copmatnew, cparmatnew))
#         }
#         return(NULL)
#     }
#     ## --- END special cases. ---
#     Ac <- GtoA(rvc$G)
#     ## Get variables that won't work:
#     bad <- unique(as.vector(Ac[1:ntrunc, (ntrunc+1):d]))
#     if (leaf %in% bad) return(NULL)
#     ## Leaf if valid. Reorder vine:
#     vc <- vars(rvc)
#     ileaf <- which(vc == leaf)
#     Ac[, c(ileaf, d)] <- Ac[, c(d, ileaf)]
#     copmatc <- rvc$copmat
#     copmatc[, c(ileaf, d)] <- copmatc[, c(d, ileaf)]
#     cparmatc <- rvc$cparmat
#     cparmatc[, c(ileaf, d)] <- cparmatc[, c(d, ileaf)]
#     rvine(AtoG(Ac), copmatc, cparmatc)
# }
