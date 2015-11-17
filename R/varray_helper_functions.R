#' Relabel Variables in a Vine Array
#'
#' @param A Vine array. Could be truncated.
#' @param labs Vector of new labels. The order of the labels correspond to
#' the order of the variables in the vine array \code{A}.
#' @details Identical to \code{CopulaModel::varrayperm} but allows for
#' the posibility that \code{A} is not square.
#' @return A relabelled vine array matrix.
#' @examples
#' (A <- truncvarray(CopulaModel::Cvinearray(5), 2))
#' relabel.varray(A, c(3, 2, 1, 5, 4))
#' @export
relabel.varray <- function(A, labs = 1:ncol(A)) {
    p <- ncol(A)
    r <- nrow(A)
    labs_orig <- varray.vars(A)
    wch_lab <- invert.perm(labs_orig)
    for (row in 1:r) {
        for (col in row:p) {
            A[row, col] <- labs[wch_lab[A[row, col]]]
        }
    }
    A
}

#' Invert a permutation
#'
#' For a permutation of a set of integers \code{1, 2, ...,p}, finds the
#' inverse permutation using the \code{\link{which}} function.
#'
#' @param perm Vector of integers in \code{{1:length(perm)}}
#' @return A vector of length \code{length(perm)} of the inverse permutation.
#' @note This function won't check whether the integers you input are
#' in the set \code{{1:length(perm)}}, but allows for length-0 entry.
#' @examples
#' perm <- c(5, 1, 2, 3, 4)
#' (perminv <- invert.perm(perm))
#' perm[perminv]
#' perminv[perm]
#'
#' ## The zero case:
#' invert.perm(integer(0))
#' @export
invert.perm <- function(perm) {
    p <- length(perm)
    if (p <= 1) return(perm)
    sapply(1:p, function(i) which(perm == i))
}

#' Truncate a Vine Array
#'
#' Truncates a vine array, collapsing the variables upwards.
#'
#' @param A A vine array, possibly truncated.
#' @param ntrunc Integer; truncation level
#' @details If \code{ntrunc >= nrow(A) - 1}, the original vine array is
#' returned. Otherwise, a truncated vine array with \code{ntrunc + 1} rows is
#' returned. The variables are listed along the initial diagonal of the vine
#' array, then continue along the bottom row.
#' @examples
#' (A <- CopulaModel::Dvinearray(6))
#' (A <- truncvarray(A, 3))
#' (A <- relabel.varray(A, c(6, 2, 4, 3, 1, 5)))
#' truncvarray(A, 2)
#' @export
truncvarray <- function(A, ntrunc) {
    p <- ncol(A)
    r <- nrow(A)
    if (ntrunc > r-2) return(A)
    vars <- varray.vars(A)
    A <- A[1:ntrunc, 1:p]
    if (!is.matrix(A)) A <- matrix(A, ncol = p)
    vars[1:ntrunc] <- 0
    rbind(A, matrix(vars, nrow=1))
}

#' Extract Variables in a Vine Array
#'
#' Extract variables in a vine array, possibly truncated, in the order
#' of the vine array.
#'
#' @param A Vine array, possibly truncated.
#' @return Vector of vine variables.
#' @examples
#' A <- CopulaModel::Dvinearray(5)
#' A <- relabel.varray(A, c(5, 2, 4, 3, 1))
#' varray.vars(A)
#'
#' A <- truncvarray(A, 2)
#' varray.vars(A)
#' @export
varray.vars <- function(A) {
    p <- ncol(A)
    r <- nrow(A)
    firstvars <- diag(A)
    secondvars <- A[r, r+seq_len(p-r)]
    c(firstvars, secondvars)
}

#' Center a Vine Array
#'
#' Converts a vine array \code{A} so that the first variables (up to
#' truncation) are not leaves. So, a slightly weaker condition than
#' natural order.
#'
#' @param A Vine array, possibly truncated.
#' @details For a \code{t}-truncated vine array \code{(t < ncol(A)-1)},
#' the vine array is re-ordered so that the first \code{t} variables
#' introduced in the outputted array are not leaves.
#'
#' If \code{t = ncol(A)-1}, then the entered vine isn't truncated, and the
#' first \code{t-1} variables in the outputted array are not leaves (in fact,
#' a natural order array is outputted, since it satisfies that requirement).
#' @export
center.varray <- function(A) {
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

#' Convert Vine Array to Convenient Array
#'
#' (A "convenient" array is achieved by moving the variables in a vine array,
#' possibly truncated, to the top row)
#' @param A Vine array
#' @param Acon A convenient vine array
#' @rdname A_Acon_convert
Atocon <- function(A) {
    ntrunc <- nrow(A) - 1
    vars <- varray.vars(A)
    Acon <- A[1:ntrunc, ]
    if (!is.matrix(Acon)) Acon <- matrix(Acon, nrow = ntrunc)
    diag(Acon) <- 0
    rbind(matrix(vars, nrow = 1), Acon)
}

#' @rdname A_Acon_convert
contoA <- function(Acon) {
    ntrunc <- nrow(Acon) - 1
    d <- ncol(Acon)
    vars <- Acon[1, ]
    A <- Acon[-1, ]
    if (!is.matrix(A)) A <- matrix(A, nrow = 1)
    A <- rbind(A, matrix(c(rep(0, ntrunc), vars[(ntrunc+1):d]), nrow = 1))
    diag(A) <- vars[1:(ntrunc+1)]
    A
}
