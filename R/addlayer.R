#' Add Layers to a Regular Vine
#'
#' @param rv Object of type "rvine".
#' @param a The vine array columns to add. Could be a vector if only one
#' variable is being added; or, could be a list of vectors for each variable
#' being added. Or it could be the matrix itself.
#' @param cops The copula models for each edge. Input as in \code{a}.
#' @param cpars The copula parameters for each copula model. Input as in \code{a}.
#' @return Object of type 'rvine' with the new layer(s) added.
#' @export
addlayer.rvine <- function(rv, a, cops, cpars) {
    if (!is.matrix(a)) a <- makevinemat(a)
    if (!is.matrix(cops)) cops <- makevinemat(cops)
    if (!is.matrix(cpars)) cpars <- makevinemat(cpars)
    ## Combine arrays:
    A <- addlayer_mat(rv$A, a, 0)
    ## Combine copmats:
    copmat <- addlayer_mat(rv$copmat, cops, "")
    ## Combine cparmats:
    cparmat <- addlayer_mat(rv$copmat, cpars, list(NULL))
    rvine(A, copmat, cparmat)
}

#' @export
addlayer <- function(...) UseMethod("addlayer")

#' Remove Upper Layers in a Vine
#'
#' @param rv Object of type "rvine"
#' @param num Integer; number of layers to remove.
#' @return Object of type "rvine".
#' @export
rmvlayer.rvine <- function(rv, num = 1) {
    var <- vars(rv)
    nvar <- length(var)
    if (num >= nvar) return(rvine(matrix(nrow=0, ncol=0)))
    if (num == nvar - 1) return(rvine(matrix(var[1])))
    ## Peel-back matrices:
    A <- rv$A[, 1:(nvar-num)]
    copmat <- rv$copmat[, 1:(nvar-num)]
    cparmat <- rv$cparmat[, 1:(nvar-num)]
    ## Get ntrunc vector
    ntrunc <- max(apply(A, 2, function(col) sum(col != 0) - 1))
    ## Remove bottom rows if necessary.
    A <- A[1:(ntrunc+1), ]
    copmat <- copmat[seq_len(ntrunc), ]
    cparmat <- cparmat[seq_len(ntrunc), ]
    if (ntrunc == 0){
        A <- matrix(A, nrow = 1)
    }
    if (ntrunc == 1){
        copmat <- matrix(copmat, nrow = 1)
        cparmat <- matrix(cparmat, nrow = 1)
    }
    rvine(A, copmat, cparmat)
}

#' @export
rmvlayer <- function(...) UseMethod("rmvlayer")

addlayer_mat <- function(basemat, layermat, zero) {
    M1 <- basemat
    M2 <- layermat
    nrow1 <- nrow(M1)
    nrow2 <- nrow(M2)
    nrow <- max(nrow1, nrow2)
    M1 <- rbind(M1, matrix(zero, ncol=ncol(M1), nrow = nrow - nrow1))
    M2 <- rbind(M2, matrix(zero, ncol=ncol(M2), nrow = nrow - nrow2))
    cbind(M1, M2)
}
