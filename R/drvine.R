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
#' the evaluated joint density of the variables in \code{G}.
#' @examples
#' set.seed(123)
#' G <- AtoG(CopulaModel::Cvinearray(4))[1:3, ]
#' copmat <- makeuppertri(c("gum", "gal", "bvtcop",
#'                          "bvncop", "frk"), row = 2, col = 4, blanks = "")
#' cparmat <- makeuppertri.list(c(1.5, 1.5, 0.9, 3, 0.1, 0.5),
#'                              len = c(1,1,2,1,1), row = 2, col = 4)
#' rv <- rvine(G=G, copmat=copmat, cparmat=cparmat)
#' dat <- rrvine(10, rv)
#' logdrvine(dat, rv)
#' drvine(runif(4), rv)
#'
#' ## The variables in G don't need to refer to all data:
#' u <- runif(4)
#' drvine(u, subset(rv, 1:2))
#' ## ...is the same as:
#' drvine(u[1:2], subset(rv, 1:2))
#' @seealso \code{\link{rvine}}
#' @rdname d_logd_rvine
#' @export
logdrvine <- function(dat, rv) {
    ## Array
    A <- GtoA(rv$G)
    nrowA <- nrow(A)
    ncolA <- ncol(A)
    if (ncolA == 0) return(NULL)
    if (nrowA == 1) return(0) # Independence vine
    ntrunc <- nrowA - 1
    copmat <- rv$copmat
    cparmat <- rv$cparmat
    if (is.vector(dat)) dat <- matrix(dat, nrow = 1)
    if (ncolA == 2) {
        logdcop <- get(paste0("logd", copmat[1, 2]))
        u <- dat[, A[1,1]]
        v <- dat[, A[2,2]]
        return(logdcop(u, v, cparmat[1, 2][[1]]))
    }
    ## Get parvec:
    if (is.list(cparmat[1,1])) {
        parvec <- c(t(cparmat), recursive = TRUE)
    } else {
        parvec <- t(cparmat)[lower.tri(t(cparmat))]
    }

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
    if (!is.matrix(dat)) dat <- matrix(dat, nrow = 1)
    ## Make array variables 1:ncol(A), and permute data to reflect that.
    v <- vars(rv)
    A <- GtoA(relabel(rv)$G) # relabel.varray(A)
    dat <- dat[, v]
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
drvine <- function(dat, rv) exp(logdrvine(dat, rv))
