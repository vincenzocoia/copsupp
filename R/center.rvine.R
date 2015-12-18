#' Center a Vine Array
#'
#' Converts a vine array \code{A} so that the first variables (up to
#' truncation) are not leaves. So, a slightly weaker condition than
#' natural order. Deprecated; use \code{\link{center.rvine}} instead.
#'
#' @param rv Regular vine object.
#' @details For a \code{t}-truncated vine array \code{(t < ncol(A)-1)},
#' the vine array is re-ordered so that the first \code{t} variables
#' introduced in the outputted array are not leaves.
#'
#' If \code{t = ncol(A)-1}, then the entered vine isn't truncated, and the
#' first \code{t-1} variables in the outputted array are not leaves (in fact,
#' a natural order array is outputted, since it satisfies that requirement).
#' @return A regular vine, with variables reordered.
#' @examples
#' ## Setup a vine.
#' A <- CopulaModel::Dvinearray(5)
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
#' (rv <- rvine(A, copmat, cparmat,
#'              list(pexp, identity, pnorm, pexp, sqrt)))
#'
#' ## Center it. Since it's a complete vine, the output is in natural order:
#' center(rv)
#'
#' ## Center a truncated vine:
#' center(trunc(rv, 2))
#'
#' ## It works in these cases too:
#' center(trunc(rv, 0))
#' center(subset(rv, 5))
#' center(subset(rv, integer(0)))
#' @export
center.rvine <- function(rv) {
    A <- rv$A
    ntrunc <- nrow(A) - 1
    ## Nothing to do if it's an independence vine:
    if (ntrunc == 0) return(rv)
    d <- ncol(A)
    ## Nothing to do if there's less than 3 variables:
    if (d == 0 | d == 2) return(rv)
    copmat <- rv$copmat
    cparmat <- rv$cparmat
    marg <- rv$marg
    ## If the vine is complete, just use "natural order":
    if (ntrunc == d-1) {
        Anew <- varray2NO(A)$NOa
        if (!is.null(copmat)) copmat <- reformcopmat(copmat, Anew = Anew, Aold = A)
        if (!is.null(cparmat)) cparmat <- reformcopmat(cparmat, Anew = Anew, Aold = A)
        vA <- vars(rv)
        vnew <- vars(rvine(Anew))
        vmap <- sapply(vA, function(v_) which(vnew == v_))
        return(rvine(Anew, copmat, cparmat, rv$marg[vmap]))
    }
    ## Initiate the final array as (ntrunc+1)x(ntrunc+1) array (call it B)
    ##  using variables in A[, d], with A[d,d] going at the end.
    Bvars <- A[, d]
    B <- subset(rv, Bvars)$A
    ovars <- vars(rv)  # Stands for "original variables"
    rvars <- setdiff(ovars, Bvars)  # Stands for "remaining variables".
    ## Convert A to a convenient form by moving labels to top row:
    Acon <- Atocon(A)
    ## Fill in B until there's no more variables left to fill:
    while (length(rvars) > 0) {
        ## Which of the remaining variables are in the next layer of the vine?
        layer <- integer(0)
        for (col in (ntrunc+1):d) {
            tf <- rvars %in% A[, col]
            if (sum(tf) == 1) layer <- c(layer, rvars[tf])
        }
        rvars <- setdiff(rvars, layer)
        wchremain <- sapply(rvars, function(i) which(ovars == i))
        if (length(wchremain) > 0) Asub <- Acon[, -wchremain] else Asub <- Acon
        if (!is.matrix(Acon)) Acon <- matrix(Acon, nrow = ntrunc + 1)
        ## Add the variables in the next layer:
        for (v in layer) {
            nextcol <- matrix(nrow = ntrunc + 1)
            nextcol[ntrunc + 1, ] <- v
            Asubsub <- Asub
            for (t in 1:ntrunc) {
                Asubsub <- Asubsub[, -1]
                if (!is.matrix(Asubsub)) Asubsub <- matrix(Asubsub, nrow = ntrunc + 1)
                ## Get the possible nodes for this tree
                nodes <- Asubsub[c(1, t+1), ]
                if (!is.matrix(nodes)) nodes <- matrix(nodes, nrow = 2)
                ## Make sure the conditioned variables are correct:
                keepcols <- rep(TRUE, ncol(Asubsub))
                if (t > 1) {
                    condn <- Asubsub[2:t, ]
                    if (!is.matrix(condn)) condn <- matrix(condn, nrow = t-1)
                    keepcols <- apply(condn, 2, function(col)
                        all(sort(col) == sort(nextcol[1:(t-1), 1])))
                }
                ## Find variable v and its partner.
                wchprsnt <- (nodes == v) & matrix(keepcols, byrow = TRUE,
                                                  ncol = ncol(Asubsub), nrow = 2)
                nextcol[t, 1] <- nodes[wchprsnt[2:1, ]]
            }
            ## Add the column to B:
            B <- cbind(B, nextcol)
        }
    }
    ## Rearrange the copula and parameter matrices, and marginals:
    if (!is.null(copmat)) copmat <- reformcopmat(copmat, B, A)
    if (!is.null(cparmat)) cparmat <- reformcopmat(cparmat, B, A)
    newvars <- vars(rvine(B))
    marg <- marg[sapply(ovars, function(v_) which(newvars == v_))]
    rvine(B, copmat, cparmat, marg)
}

#' @export
center <- function(...) UseMethod("center")

#' Center a Vine Array
#'
#' Converts a vine array \code{A} so that the first variables (up to
#' truncation) are not leaves. So, a slightly weaker condition than
#' natural order. Deprecated; use \code{\link{center.rvine}} instead.
#'
#' @param A Vine array, possibly truncated.
#' @details For a \code{t}-truncated vine array \code{(t < ncol(A)-1)},
#' the vine array is re-ordered so that the first \code{t} variables
#' introduced in the outputted array are not leaves.
#'
#' If \code{t = ncol(A)-1}, then the entered vine isn't truncated, and the
#' first \code{t-1} variables in the outputted array are not leaves (in fact,
#' a natural order array is outputted, since it satisfies that requirement).
#' @return A vine array with the same dimensions as \code{A}.
#' @export
centervarray <- function(A) {
    warning("'centervarray()' is deprecated. Please use 'center.rvine()' instead.")
    ntrunc <- nrow(A) - 1
    d <- ncol(A)
    if (ntrunc == d-1) return(CopulaModel::varray2NO(A)$NOa)
    ## Initiate the final array as (ntrunc+1)x(ntrunc+1) array (call it B)
    ##  using variables in A[, d], with A[d,d] going at the end.
    Bvars <- A[, d]
    B <- rvinesubset(A, Bvars)
    ovars <- varray.vars(A)  # Stands for "original variables"
    rvars <- setdiff(ovars, Bvars)  # Stands for "remaining variables".
    ## Convert A to a convenient form by moving labels to top row:
    Acon <- Atocon(A)
    ## Fill in B until there's no more variables left to fill:
    while (length(rvars) > 0) {
        ## Which of the remaining variables are in the next layer of the vine?
        layer <- integer(0)
        for (col in (ntrunc+1):d) {
            tf <- rvars %in% A[, col]
            if (sum(tf) == 1) layer <- c(layer, rvars[tf])
        }
        rvars <- setdiff(rvars, layer)
        wchremain <- sapply(rvars, function(i) which(ovars == i))
        if (length(wchremain) > 0) Asub <- Acon[, -wchremain] else Asub <- Acon
        if (!is.matrix(Acon)) Acon <- matrix(Acon, nrow = ntrunc + 1)
        ## Add the variables in the next layer:
        for (v in layer) {
            nextcol <- matrix(nrow = ntrunc + 1)
            nextcol[ntrunc + 1, ] <- v
            Asubsub <- Asub
            for (t in 1:ntrunc) {
                Asubsub <- Asubsub[, -1]
                if (!is.matrix(Asubsub)) Asubsub <- matrix(Asubsub, nrow = ntrunc + 1)
                ## Get the possible nodes for this tree
                nodes <- Asubsub[c(1, t+1), ]
                if (!is.matrix(nodes)) nodes <- matrix(nodes, nrow = 2)
                ## Make sure the conditioned variables are correct:
                keepcols <- rep(TRUE, ncol(Asubsub))
                if (t > 1) {
                    condn <- Asubsub[2:t, ]
                    if (!is.matrix(condn)) condn <- matrix(condn, nrow = t-1)
                    keepcols <- apply(condn, 2, function(col)
                        all(sort(col) == sort(nextcol[1:(t-1), 1])))
                }
                ## Find variable v and its partner.
                wchprsnt <- (nodes == v) & matrix(keepcols, byrow = TRUE,
                                                  ncol = ncol(Asubsub), nrow = 2)
                nextcol[t, 1] <- nodes[wchprsnt[2:1, ]]
            }
            ## Add the column to B:
            B <- cbind(B, nextcol)
        }
    }
    B
}
