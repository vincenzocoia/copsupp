#' Center a Vine
#'
#' Converts a vine array (\code{centervarray}) or vine
#' (\code{center.rvine}) so that the first variables (up to
#' truncation) are not "leaves". This may not work outside of
#' regular vines with truncation under the traditional sense.
#'
#' @param rv Regular vine object.
#' @param G Vine array. Must be a regular vine, with truncation in
#' the traditional sense.
#' @details For a \code{t}-truncated vine array \code{(t < ncol(G)-1)},
#' the vine array is re-ordered so that the first \code{t} variables
#' introduced in the outputted array are not leaves.
#'
#' If \code{t = ncol(G)-1}, then the entered vine isn't truncated, and the
#' first \code{t-1} variables in the outputted array are not leaves (in fact,
#' a natural order array is outputted, since it satisfies that requirement).
#' @return G-vine array or vine, with variables reordered.
#' @examples
#' ## Setup a vine.
#' G <- AtoG(CopulaModel::Dvinearray(5))
#' copmat <- makeuppertri(c("gum", "mtcj", "gal", "joe",
#'                          "frk", "gum", "bb7",
#'                          "bb1", "indepcop",
#'                          "bb8"), 4, 5, "")
#' cparmat <- makeuppertri.list(c(3, 2.5, 2, 1.5,
#'                                1, 1.3, 2, 2,
#'                                3, 4,
#'                                5, 0.5),
#'                              len = c(1,1,1,1,1,1,2,2,0,2),
#'                              4, 5)
#' (rv <- rvine(G, copmat, cparmat))
#'
#' ## Center it. Since it's a complete vine, the output is in natural order:
#' centervarray(G)
#' center(rv)
#'
#' ## Center a truncated vine:
#' centervarray(G[1:3, ])
#' center(trunc(rv, 2))
#'
#' ## It works in these cases too:
#' center(trunc(rv, 0))
#' center(subset(rv, 5))
#' center(subset(rv, integer(0)))
#' @import CopulaModel
#' @rdname center
#' @export
centervarray <- function(G) {
    ntrunc <- nrow(G) - 1
    ## Nothing to do if it's an independence vine:
    if (ntrunc == 0) return(G)  # Accounts for one-variable case too.
    d <- ncol(G)
    ## Nothing to do if there's less than 3 variables:
    if (d <= 2) return(G)
    ovars <- G[1, ]
    ## If the vine is complete, just use "natural order":
    if (ntrunc == d-1) {
        ## Change variable names to variable order, and convert to natural order.
        Gnew <- relabelvarray(G)
        Gnew <- AtoG(varray2NO(GtoA(Gnew))$NOa)
        ## Convert back to variable names:
        vnew <- ovars[Gnew[1, ]]
        Gnew <- relabelvarray(Gnew, vnew)
        return(Gnew)
    }
    ## Initiate a "centralized" order of the variables.
    ord <- G[c(2:(ntrunc+1), 1), d]
    ## Candidates for the next variables:
    candvars <- setdiff(G[1, ], ord)
    while (length(candvars) > 0) {
        ## Collect variables to go into the next "shell" of the vine.
        nextshell <- integer(0)
        for (j in (ntrunc+1):d) {
            x <- intersect(G[, j], candvars)
            if (length(x) == 1) nextshell <- c(nextshell, x)
        }
        ## Remove the new shell variables from the candidates, and add them
        ##  to the order
        ord <- c(ord, nextshell)
        candvars <- setdiff(candvars, nextshell)
    }
    ## Get re-ordered vine array
    reordervarray(G, ord)
}

#' @rdname center
#' @export
center.rvine <- function(rv) {
    ## Center vine array first.
    G <- rv$G
    Gnew <- centervarray(G)
    ## Reform copmat and cparmat:
    if (nrow(G) > 1) {
        copmat <- reformcopmat(rv$copmat, Gnew, G)
        cparmat <- reformcopmat(rv$cparmat, Gnew, G)
        return(rvine(Gnew, copmat, cparmat))
    } else {
        return(rvine(Gnew))
    }
}

#' @export
center <- function(...) UseMethod("center")


#
#     ## Initiate the final array as (ntrunc+1)x(ntrunc+1) array (call it B)
#     ##  using variables in G[, d], with G[1,d] going at the end.
#     Bvars <- G[, d]
#     B <- subsetvarray(G, Bvars)
#     rvars <- setdiff(ovars, Bvars)  # Stands for "remaining variables".
#     ## Fill in B until there's no more variables left to fill:
#     while (length(rvars) > 0) {
#         ## Which of the remaining variables are in the next layer of the vine?
#         layer <- integer(0)
#         for (col in (ntrunc+1):d) {
#             tf <- rvars %in% G[, col]
#             if (sum(tf) == 1) layer <- c(layer, rvars[tf])
#         }
#         rvars <- setdiff(rvars, layer)
#         wchremain <- sapply(rvars, function(i) which(ovars == i))
#         if (length(wchremain) > 0) Asub <- A[, -wchremain] else Asub <- A
#         if (!is.matrix(A)) A <- matrix(A, nrow = ntrunc + 1)
#         ## Add the variables in the next layer:
#         for (v in layer) {
#             nextcol <- matrix(nrow = ntrunc + 1)
#             nextcol[ntrunc + 1, ] <- v
#             Asubsub <- Asub
#             for (t in 1:ntrunc) {
#                 Asubsub <- Asubsub[, -1]
#                 if (!is.matrix(Asubsub)) Asubsub <- matrix(Asubsub, nrow = ntrunc + 1)
#                 ## Get the possible nodes for this tree
#                 nodes <- Asubsub[c(1, t+1), ]
#                 if (!is.matrix(nodes)) nodes <- matrix(nodes, nrow = 2)
#                 ## Make sure the conditioned variables are correct:
#                 keepcols <- rep(TRUE, ncol(Asubsub))
#                 if (t > 1) {
#                     condn <- Asubsub[2:t, ]
#                     if (!is.matrix(condn)) condn <- matrix(condn, nrow = t-1)
#                     keepcols <- apply(condn, 2, function(col)
#                         all(sort(col) == sort(nextcol[1:(t-1), 1])))
#                 }
#                 ## Find variable v and its partner.
#                 wchprsnt <- (nodes == v) & matrix(keepcols, byrow = TRUE,
#                                                   ncol = ncol(Asubsub), nrow = 2)
#                 nextcol[t, 1] <- nodes[wchprsnt[2:1, ]]
#             }
#             ## Add the column to B:
#             B <- cbind(B, nextcol)
#         }
#     }
#     ## Rearrange the copula and parameter matrices, and marginals:
#     AtoG(B)
# }
#
#
#
# center.rvine <- function(rv) {
#     G <- rv$G
#     ntrunc <- nrow(G) - 1
#     ## Nothing to do if it's an independence vine:
#     if (ntrunc == 0) return(rv)
#     d <- ncol(G)
#     ## Nothing to do if there's less than 3 variables:
#     if (d == 0 | d == 2) return(rv)
#     copmat <- rv$copmat
#     cparmat <- rv$cparmat
#     ovars <- vars(rv)
#     ## If the vine is complete, just use "natural order":
#     if (ntrunc == d-1) {
#         ## Change variable names to variable order, and convert to natural order.
#         ord2var <- vars(rv)
#         rv <- relabel(rv, 1:d)
#         Gnew <- AtoG(varray2NO(GtoA(G))$NOa)
#         ## Convert back to variable names:
#         newordvars <- Gnew[1, ]
#         vnew <- ord2var[newordvars]
#         Gnew <- relabelvarray(Gnew, vnew)
#         if (!is.null(copmat)) copmat <- reformcopmat(copmat, Gnew = Gnew, Gold = G)
#         if (!is.null(cparmat)) cparmat <- reformcopmat(cparmat, Gnew = Gnew, Gold = G)
#         vmap <- sapply(vnew, function(vnew_) which(ovars == vnew_))
#         return(rvine(Gnew, copmat, cparmat, rv$marg[vmap]))
#     }
#     ## Initiate the final array as (ntrunc+1)x(ntrunc+1) array (call it B)
#     ##  using variables in G[, d], with G[d,d] going at the end.
#     Bvars <- G[, d]
#     B <- subset(rv, Bvars)$G
#     rvars <- setdiff(ovars, Bvars)  # Stands for "remaining variables".
#     ## Fill in B until there's no more variables left to fill:
#     while (length(rvars) > 0) {
#         ## Which of the remaining variables are in the next layer of the vine?
#         layer <- integer(0)
#         for (col in (ntrunc+1):d) {
#             tf <- rvars %in% G[, col]
#             if (sum(tf) == 1) layer <- c(layer, rvars[tf])
#         }
#         rvars <- setdiff(rvars, layer)
#         wchremain <- sapply(rvars, function(i) which(ovars == i))
#         if (length(wchremain) > 0) Gsub <- G[, -wchremain] else Gsub <- G
#         if (!is.matrix(G)) G <- matrix(G, nrow = ntrunc + 1)
#         ## Add the variables in the next layer:
#         for (v in layer) {
#             nextcol <- matrix(nrow = ntrunc + 1)
#             nextcol[ntrunc + 1, ] <- v
#             Gsubsub <- Gsub
#             for (t in 1:ntrunc) {
#                 Gsubsub <- Gsubsub[, -1]
#                 if (!is.matrix(Gsubsub)) Gsubsub <- matrix(Gsubsub, nrow = ntrunc + 1)
#                 ## Get the possible nodes for this tree
#                 nodes <- Gsubsub[c(1, t+1), ]
#                 if (!is.matrix(nodes)) nodes <- matrix(nodes, nrow = 2)
#                 ## Make sure the conditioned variables are correct:
#                 keepcols <- rep(TRUE, ncol(Gsubsub))
#                 if (t > 1) {
#                     condn <- Gsubsub[2:t, ]
#                     if (!is.matrix(condn)) condn <- matrix(condn, nrow = t-1)
#                     keepcols <- apply(condn, 2, function(col)
#                         all(sort(col) == sort(nextcol[1:(t-1), 1])))
#                 }
#                 ## Find variable v and its partner.
#                 wchprsnt <- (nodes == v) & matrix(keepcols, byrow = TRUE,
#                                                   ncol = ncol(Gsubsub), nrow = 2)
#                 nextcol[t, 1] <- nodes[wchprsnt[2:1, ]]
#             }
#             ## Add the column to B:
#             B <- cbind(B, nextcol)
#         }
#     }
#     ## Rearrange the copula and parameter matrices, and marginals:
#     if (!is.null(copmat)) copmat <- reformcopmat(copmat, B, G)
#     if (!is.null(cparmat)) cparmat <- reformcopmat(cparmat, B, G)
#     newvars <- vars(rvine(B))
#     rvine(B, copmat, cparmat)
# }
#
#
#
# #' Center a Vine Array
# #'
# #' Converts a vine array \code{A} so that the first variables (up to
# #' truncation) are not leaves. So, a slightly weaker condition than
# #' natural order. Deprecated; use \code{\link{center.rvine}} instead.
# #'
# #' @param A Vine array, possibly truncated.
# #' @details For a \code{t}-truncated vine array \code{(t < ncol(A)-1)},
# #' the vine array is re-ordered so that the first \code{t} variables
# #' introduced in the outputted array are not leaves.
# #'
# #' If \code{t = ncol(A)-1}, then the entered vine isn't truncated, and the
# #' first \code{t-1} variables in the outputted array are not leaves (in fact,
# #' a natural order array is outputted, since it satisfies that requirement).
# #' @return A vine array with the same dimensions as \code{A}.
# #' @export
# centervarray <- function(A) {
#     warning("'centervarray()' is deprecated. Please use 'center.rvine()' instead.")
#     ntrunc <- nrow(A) - 1
#     d <- ncol(A)
#     if (ntrunc == d-1) return(CopulaModel::varray2NO(A)$NOa)
#     ## Initiate the final array as (ntrunc+1)x(ntrunc+1) array (call it B)
#     ##  using variables in A[, d], with A[d,d] going at the end.
#     Bvars <- A[, d]
#     B <- rvinesubset(A, Bvars)
#     ovars <- varray.vars(A)  # Stands for "original variables"
#     rvars <- setdiff(ovars, Bvars)  # Stands for "remaining variables".
#     ## Convert A to a convenient form by moving labels to top row:
#     Acon <- Atocon(A)
#     ## Fill in B until there's no more variables left to fill:
#     while (length(rvars) > 0) {
#         ## Which of the remaining variables are in the next layer of the vine?
#         layer <- integer(0)
#         for (col in (ntrunc+1):d) {
#             tf <- rvars %in% A[, col]
#             if (sum(tf) == 1) layer <- c(layer, rvars[tf])
#         }
#         rvars <- setdiff(rvars, layer)
#         wchremain <- sapply(rvars, function(i) which(ovars == i))
#         if (length(wchremain) > 0) Asub <- Acon[, -wchremain] else Asub <- Acon
#         if (!is.matrix(Acon)) Acon <- matrix(Acon, nrow = ntrunc + 1)
#         ## Add the variables in the next layer:
#         for (v in layer) {
#             nextcol <- matrix(nrow = ntrunc + 1)
#             nextcol[ntrunc + 1, ] <- v
#             Asubsub <- Asub
#             for (t in 1:ntrunc) {
#                 Asubsub <- Asubsub[, -1]
#                 if (!is.matrix(Asubsub)) Asubsub <- matrix(Asubsub, nrow = ntrunc + 1)
#                 ## Get the possible nodes for this tree
#                 nodes <- Asubsub[c(1, t+1), ]
#                 if (!is.matrix(nodes)) nodes <- matrix(nodes, nrow = 2)
#                 ## Make sure the conditioned variables are correct:
#                 keepcols <- rep(TRUE, ncol(Asubsub))
#                 if (t > 1) {
#                     condn <- Asubsub[2:t, ]
#                     if (!is.matrix(condn)) condn <- matrix(condn, nrow = t-1)
#                     keepcols <- apply(condn, 2, function(col)
#                         all(sort(col) == sort(nextcol[1:(t-1), 1])))
#                 }
#                 ## Find variable v and its partner.
#                 wchprsnt <- (nodes == v) & matrix(keepcols, byrow = TRUE,
#                                                   ncol = ncol(Asubsub), nrow = 2)
#                 nextcol[t, 1] <- nodes[wchprsnt[2:1, ]]
#             }
#             ## Add the column to B:
#             B <- cbind(B, nextcol)
#         }
#     }
#     B
# }
