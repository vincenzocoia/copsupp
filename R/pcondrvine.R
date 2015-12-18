#' Conditional Distribution in a Regular Vine
#'
#' Evaluates the conditional distribution of a variable in a regular vine,
#' given values of the other variables in that vine.
#'
#' @param dat vector or matrix of observations (columns are variables).
#' @param rv Regular vine object, complete.
#' @param cond Integer; the variable you wish to condition on (i.e. the
#' column number of \code{dat}, also present in \code{rv}).
#' @param vbls Vector of integers; the variables you wish to consider. Default
#' is all variables in \code{rv}.
#' @param verbose Logical; should the function output how it goes about
#' finding the conditional distribution?
#' @details To compute the conditional distribution, the vine is subsetted to
#' the selected variables if possible. Then, if the conditioned variable is a
#' leaf, the conditional distributon is directly computed. If it's not a leaf,
#' the conditional distribution is computed by integrating the density.
#'
#' If the subsetted vine does not exist, then the vine will be subsetted "as
#' much as possible", and the remaining variables that cannot be removed
#' will be integrated out to find the joint density of the selected variables,
#' from which the conditional cdf will be found.
#' @return A vector of length = the number of observations in \code{dat},
#' representing the evaluated conditional distribution of variable \code{cond}
#' given the other variables in \code{vbls}.
#' @examples
#' ## D-Vine example
#' A <- CopulaModel::Dvinearray(5)
#' A <- relabel.varray(A, c(1, 5, 4, 3, 2))
#' A <- truncvarray(A, 2)
#' copmat <- makeuppertri("bvncop", 2, 5)
#' cparmat <- makeuppertri(c(1:7/10), 2, 5, byRow = FALSE)
#' udat <- fvinesim(10, A, copmat, cparmat)
#' pcondrvine(udat, 5, A, copmat, cparmat, verbose=T)  # integrates vine density.
#' pcondrvine(udat, 2, A, copmat, cparmat, verbose=T)  # computes from D-vine formula
#'
#' ## C-Vine example
#' A <- CopulaModel::Cvinearray(5)
#' A <- truncvarray(A, 2)
#' udat <- fvinesim(10, A, copmat, cparmat)
#' pcondrvine(udat, 3, A, copmat, cparmat, verbose=T)  # computes from general R-vine algo
#'
#' ## Array doesn't have to involve all data:
#' A <- CopulaModel::Dvinearray(5)
#' A <- truncvarray(A, 2)
#' A <- rvinesubset(A, 3:5)
#' copmat <- makeuppertri("frk", 2, 3, "")
#' cparmat <- makeuppertri(3:1, 2, 3)
#' pcondrvine(1:5/10, 3, A, copmat, cparmat)
#' pcondrvine(1:5/10, 4, A, copmat, cparmat)
#' ## are the same as...
#' A <- CopulaModel::Dvinearray(3)
#' A <- truncvarray(A, 2)
#' pcondrvine(3:5/10, 1, A, copmat, cparmat)
#' pcondrvine(3:5/10, 2, A, copmat, cparmat)
#' @export
pcondrvine <- function(dat, rv, cond, vbls = vars(rv), verbose = FALSE) {
    if (is.vector(dat)) dat <- matrix(dat, nrow = 1)
    ## Extract info
    A <- rv$A
    Fmarg <- rv$marg
    d <- ncol(A)
    ptot <- ncol(dat)
    ntrunc <- nrow(A) - 1
    v <- vars(rv)
    ikeep <- sapply(vbls, function(vbl) which(v == vbl))
    ## Uniformize data:
    dat <- dat[, vbls]
    for (i in 1:length(ikeep)) {
        dat[, i] <- Fmarg[[ikeep[i]]](dat[, i])
    }
    ## Can I subset the vine?
    subrv <- subset(rv, vbls)
    if (!is.null(subrv)) {

    }
    ## Is cond a leaf? If so, get the vine array with it as a leaf.
    Aleaf <- releaf.varray(A, leaf = cond)
    if (is.null(Aleaf)) {
        if (verbose) print(paste0("cond=", cond, " is not a leaf. ",
                                 "Obtaining conditional distribution by integration."))
        res <- apply(udat, 1, function(row) {
            dens <- function(xcond){  # Accepts uniform variable 'cond'.
                x <- row
                x[cond] <- xcond
                drvine(x, A, copmat, cparmat)
            }
            integrate.mv(dens, 0, row[cond]) / integrate.mv(dens, 0, 1)
        })
    } else {
        ## We'll need to re-arrange the copmat and cparmat to match Aleaf.
        copmat <- reform.copmat(copmat, Aleaf, A)
        cparmat <- reform.copmat(cparmat, Aleaf, A)
        ## We'll need to re-arrange the data so that it's in order of the
        ##  vine array.
        revars <- varray.vars(Aleaf)
        udat <- udat[, revars]
        if (!is.matrix(udat)) udat <- matrix(udat, ncol = d)
        if (ncol(udat) <= 2) {
            ## NOTE: This special case is needed because
            ##  rVineTruncCondCDF() won't accept a vine with 2 variables.
            if (verbose) print(paste0("cond=", cond, " is one of a pair. ",
                                      "Using `pcondcop()`."))
            library(CopulaModel)
            pcondcop <- get(paste0("pcond", copmat[1, 2]))
            cpar <- cparmat[1, 2]
            if (is.list(cpar)) cpar <- cpar[[1]]
            res <- apply(udat, 1, function(row) pcondcop(row[2], row[1], cpar))
        } else {
            if (verbose) print(paste0("cond=", cond, " is a leaf. ",
                                     "Using `copreg::rVineTruncCondCDF()`."))
            ## Use Bo's function
            ## Re-label the variables in the vine so that they're 1:d, in that
            ##  order.
            Aleaf <- relabel.varray(Aleaf)
            ## Fill-in vine array so it's d x d
            Aleaf <- rbind(Aleaf, matrix(0, nrow = d - ntrunc - 1, ncol = d))
            diag(Aleaf) <- 1:d
            ## parvec
            # parvec <- c(t(cparmat), recursive = TRUE)
            if (is.list(cparmat[1,1])) {
                parvec <- c(t(cparmat), recursive = TRUE)
            } else {
                parvec <- t(cparmat)[lower.tri(t(cparmat))]
            }
            ## pcondmat
            pcondmat <- apply(copmat, 1:2, function(cop) paste0("pcond", cop))
            pcondmat[!upper.tri(pcondmat)] <- ""
            ## np
            if (is.list(cparmat[1,1])) {
                np <- apply(cparmat, 1:2, function(cpar) length(cpar[[1]]))
            } else {
                np <- makeuppertri(1, nrow(cparmat), ncol(cparmat))
            }
            ## Call:
            library(CopulaModel)
            res <- copreg::rVineTruncCondCDF(parvec = parvec,
                                             udat = udat,
                                             A = Aleaf,
                                             ntrunc = ntrunc,
                                             pcondmat = pcondmat,
                                             np = np)
       }
    }
    res
}
