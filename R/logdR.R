#' Density of a Regular Vine
#'
#' Evaluates the density of a regular vine model (\code{dR}) or log density
#' (\code{logdR}).
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
#' @details This function is a wrapper for
#' \code{rvinellkv.trunc2} in the \code{CopulaModel} package.
#' @examples
#' set.seed(123)
#' A <- truncvarray(CopulaModel::Cvinearray(4), 2)
#' copmat <- makeuppertri(c("gum", "gal", "bvtcop",
#'                          "bvncop", "frk"), row = 2, col = 4, blanks = "")
#' cparmat <- makeuppertri.list(c(1.5, 1.5, 0.9, 3, 0.1, 0.5),
#'                              len = c(1,1,2,1,1), row = 2, col = 4)
#' dat <- fvinesim(10, A, copmat, cparmat)
#' logdR(dat, A, copmat, cparmat)
#' dR(c(0.5,0.5,0.5,0.5), A, copmat, cparmat)
#'
#' ## The variables in A don't need to refer to all data:
#' A <- CopulaModel::Dvinearray(6)
#' A <- rvinesubset(A, 3:6)
#' copmat <- makeuppertri("frk", 4, 4, "")
#' cparmat <- makeuppertri(6:1, 4, 4)
#' logdR(1:6/10, A, copmat, cparmat)
#' ## is the same as...
#' A <- CopulaModel::Dvinearray(4)
#' logdR(3:6/10, A, copmat, cparmat)
#' @rdname d_logd_rvine
#' @export
logdR <- function(dat, A, copmat, cparmat, Fmarg = identity) {
    ## Get parvec:
    if (is.list(cparmat[1,1])) {
        parvec <- c(t(cparmat), recursive = TRUE)
    } else {
        parvec <- t(cparmat)[lower.tri(t(cparmat))]
    }
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
        np <- matrix(0, nrow = nrow(cparmat), ncol = ncol(cparmat))
        np[upper.tri(np)] <- 1
    }
    ## Get udat
    if (length(Fmarg) == 1) Fmarg <- rep(list(Fmarg), ncolA)
    if (!is.matrix(dat)) dat <- matrix(dat, nrow = 1)
    for (col in 1:ncolA) dat[, col] <- Fmarg[[col]](dat[, col])
    ## Make array variables 1:ncol(A), and permute data to reflect that.
    vars <- varray.vars(A)
    A <- relabel.varray(A)
    dat <- dat[, vars]
    if (!is.matrix(dat)) dat <- matrix(dat, nrow = 1)
    ## Fill in A:
    A[nrowA, ] <- 0
    A <- rbind(A, matrix(0, nrow = ncolA - nrowA, ncol = ncolA))
    diag(A) <- 1:ncol(A)
    ## Use CopulaModel function
    CopulaModel::rvinellkv.trunc2(parvec, dat, A, ntrunc,
                                  logdcopmat = logdcopmat,
                                  pcondmat = pcondmat,
                                  np = np)
}

#' @rdname d_logd_rvine
#' @export
dR <- function(dat, A, copmat, cparmat, Fmarg = identity)
    exp(logdR(dat, A, copmat, cparmat, Fmarg = identity))
