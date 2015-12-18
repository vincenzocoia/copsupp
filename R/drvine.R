#' Density of a Regular Vine
#'
#' Evaluates the density of a regular vine model (\code{drvine}) or log density
#' (\code{logdrvine}).
#'
#' @param dat Data matrix. Rows are observations, and columns are variables.
#' Could be a vector if there's only one observation.
#' @param rv Regular vine object
#' @details This function is a wrapper for
#' \code{rvinellkv.trunc2} in the \code{CopulaModel} package.
#' @return Vector of length = number of observations in \code{dat}, representing
#' the evaluated joint density of the variables in \code{A}.
#' @examples
#' set.seed(123)
#' A <- truncvarray(CopulaModel::Cvinearray(4), 2)
#' copmat <- makeuppertri(c("gum", "gal", "bvtcop",
#'                          "bvncop", "frk"), row = 2, col = 4, blanks = "")
#' cparmat <- makeuppertri.list(c(1.5, 1.5, 0.9, 3, 0.1, 0.5),
#'                              len = c(1,1,2,1,1), row = 2, col = 4)
#' rv <- rvine(A=A, copmat=copmat, cparmat=cparmat)
#' dat <- fvinesim(10, rv)
#' logdrvine(dat, rv)
#' drvine(runif(4), rv)
#'
#' ## The variables in A don't need to refer to all data:
#' u <- runif(4)
#' drvine(u, subset(rv, 3:4))
#' ## ...is the same as:
#' drvine(u[3:4]), subset(rv, 3:4)
#' @seealso \code{\link{rvine}}
#' @rdname d_logd_rvine
#' @export
logdrvine <- function(dat, rv) {
    A <- rv$A
    copmat <- rv$copmat
    cparmat <- rv$cparmat
    marg <- rv$marg
    if (is.null(copmat) | is.null(cparmat))
        stop("rvine must have copmat and cparmat specified")
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
    if (is.null(marg)) marg <- rep(list(identity), ncolA)
    if (length(marg) == 1) marg <- rep(list(marg), ncolA)
    if (!is.matrix(dat)) dat <- matrix(dat, nrow = 1)
    for (col in 1:ncolA) dat[, col] <- marg[[col]](dat[, col])
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
drvine <- function(dat, rv)
    exp(logdrvine(dat, rv))
