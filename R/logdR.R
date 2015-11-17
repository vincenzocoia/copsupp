#' Density of a Regular Vine
#'
#' @param dat Data matrix. Rows are observations, and columns are variables.
#' Could be a vector if there's only one observation.
#' @param A Vine array. Integer labels should correspond to the column number
#' in \code{dat}.
#' @param copmat Matrix of bivariate copula model names, corresponding to
#' the upper-triangular part of \code{A}.
#' @param cparmat Matrix of parameters for copula models in \code{copmat}.
#' If a copula model has more or less parameters than 1, put the vector of
#' copula parameters in a list in each entry.
#' @param Fmarg List of vectorized marginal cdfs of data, in the order listed
#' in \code{dat}. Or a single such function if they're all the same.
logdR <- function(dat, A, copmat, cparmat, Fmarg = identity) {
    ## Get parvec:
    parvec <- c(cparmat, recursive = TRUE)
    ## ntrunc:
    nrowA <- nrow(A)
    ncolA <- ncol(A)
    ntrunc <- nrowA - 1
    ## logdcopmat:
    logdcopmat <- apply(copmat, 1:2, function(cop) paste0("logd", cop))
    logdcopmat[!upper.tri(logdcopmat)] <- ""
    ## pcondmat:
    pcondmat <- apply(copmat, 1:2, function(cop) paste0("pcond", cop))
    pcondmat[!upper.tri(pcondmat)] <- ""
    ## np:
    if (is.list(cparmat[1,1])) {
        np <- apply(cparmat, 1:2, function(t) length(t[[1]]))
    } else {
        np <- apply(cparmat, 1:2, length)
    }
    rvinellkv.trunc2()
    ## Fill in A:
    vars <- varray.vars(A)
    A[nrowA, ] <- 0
    A <- rbind(A, matrix(0, nrow = ncolA - nrowA, ncol = ncolA))
    diag(A) <- vars
    ## Get udat
    if (length(Fmarg) == 1) Fmarg <- rep(list(Fmarg), ncolA)
    if (!is.matrix(dat)) dat <- matrix(dat, nrow = 1)
    for (col in 1:ncolA) dat[, col] <- Fmarg[[col]](dat[, col])
    ## Use CopulaModel function
    CopulaModel::rvinellkv.trunc2(parvec, dat, A, ntrunc,
                                  logdcopmat = logdcopmat,
                                  pcondmat = pcondmat,
                                  np = np)
}
