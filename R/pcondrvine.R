#' Conditional Distribution in a Regular Vine
#'
#' Evaluates the conditional distribution of a variable in a regular vine
#' out of all specified variables. As of now, will likely only work if the vine
#' is regular, and truncated in the traditional sense.
#'
#' @param dat vector or matrix of observations (columns are variables).
#' @param rv Regular vine object.
#' @param var Integer; the variable you wish to condition on (i.e. the
#' column number of \code{dat}, also present in \code{rv}).
#' @param condset Vector of the conditioning variables (integers). Optional;
#' leave blank to condition \code{var} on all other variables in \code{rv}.
#' @param maxint Integer; maximum dimension of integration to tolerate. Put
#' \code{Inf} if you don't want an upper limit.
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
#' representing the evaluated conditional distribution of variable \code{var}
#' given the other variables in \code{condset}.
#'
#' If integration beyond \code{maxint} dimensions is required to obtain
#' the quantities, then an error is thrown.
#' @examples
#' ## D-Vine example
#' G <- AtoG(CopulaModel::Dvinearray(5))[1:3, ]
#' rv <- rvine(G, "bvncop", makeuppertri(c(1:7/10), 2, 5, byRow = FALSE))
#' dat <- rrvine(10, rv)
#'
#' ## Compute 5|1:4. There's an algorithm for that.
#' pcondrvine(dat, rv, var=5, verbose=T)
#'
#' ## Compute 5|4. pcondrvine just uses 'pcondcop()'.
#' pcondrvine(dat, rv, var=5, condset=4, verbose=T)
#'
#' ## Compute 5|2:3. Two integrals takes ~13 min if maxint > 1.
#' pcondrvine(dat, rv, var=5, condset=c(2,3), maxint=1, verbose=T)
#'
#' ## Compute 4|(1,2,3,5). No algorithm for that.
#' pcondrvine(dat, rv, var=4, verbose=T)
#' pcondrvine(dat, rv, var=4, maxint=0, verbose=T) # No int. tolerance
#' @export
pcondrvine <- function(dat, rv, var, condset, maxint = 2, verbose = FALSE) {
    v <- vars(rv)
    if (missing(condset) || is.null(condset)) {
        condset <- setdiff(v, var)
    }
    if (var %in% condset)
        stop("Conditioned variable 'var' cannot be part of conditioning variables 'condset'.")
    if (is.vector(dat)) dat <- matrix(dat, nrow = 1)
    ## Extract info
    A <- GtoA(rv$G)
    d <- length(condset) + 1
    vbls <- c(condset, var)
    ptot <- ncol(dat)
    ntrunc <- nrow(A) - 1
    ikeep <- sapply(vbls, function(vbl) which(v == vbl))
    ## Case 0: vbls = var. Just return the cdf of var.
    if (d == 1) {
        return(dat[, var])
    }
    ## Can I subset the vine?
    subrv <- subset(rv, vbls)
    if (!is.null(subrv)) {
        ## Case 1: Vine can be subsetted to 'vbls'.
        if (verbose) cat(paste0("Vine subsetted to requested variables ",
                                 paste(vbls, collapse = ", "), "\n"))
        ## Is var a leaf? If so, get the vine array with it as a leaf.
        subrvleaf <- releaf(subrv, leaf = var)
        if (is.null(subrvleaf)) {
            ## Case 1a: Subsetted vine doesn't have 'var' as a leaf.
            if (verbose) cat(paste0("var=", var, " is not a leaf."))
            if (maxint == 0)
                stop(paste("Computing 'pcondrvine' requires 1 integral.",
                            "Raise integration tolerance through 'maxint' argument."))
            if (verbose) cat("Obtaining conditional distribution by integration.\n")
            res <- apply(dat, 1, function(row) {
                dens <- function(xcond){  # Accepts uniform variable 'var'.
                    x <- row
                    x[var] <- xcond
                    drvine(x, subrv)
                }
                integrate.mv(dens, 0, row[var]) / integrate.mv(dens, 0, 1)
            })
            return(res)
        } else {
            ## Case 1b: Subsetted vine does have 'var' as a leaf.
            Aleaf <- GtoA(subrvleaf$G)
            copmat <- subrvleaf$copmat
            cparmat <- subrvleaf$cparmat
            ## We'll need to re-arrange the data so that it's in order of the
            ##  vine array.
            revars <- vars(subrvleaf)
            dat <- dat[, revars]
            subrvleaf <- relabel(subrvleaf, 1:d)
            if (ncol(dat) == 2) {
                ## Case 1ba: There's only two variables in the vine.
                ## (rVineTruncCondCDF() won't accept a vine with 2 variables)
                if (verbose) cat(paste0("var=", var, " is one of a pair. ",
                                          "Using `pcondcop()`.\n"))
                pcondcop <- get(paste0("pcond", copmat[1, 2]))
                cpar <- cparmat[1, 2][[1]]
                return(apply(dat, 1, function(row) pcondcop(row[2], row[1], cpar)))
            } else {
                ## Case 1bb: There are more than 2 variables in the vine.
                if (verbose) cat(paste0("var=", var, " is a leaf. ",
                                          "Using `copreg::rVineTruncCondCDF()`.\n"))
                ## Use Bo's function
                Aleaf <- GtoA(subrvleaf$G)
                ## Fill-in vine array so it's d x d
                Aleaf <- rbind(Aleaf, matrix(0, nrow = d - nrow(Aleaf), ncol = d))
                diag(Aleaf) <- 1:d
                ## Has problems with the independence copula: trick the function
                ##  by putting in a Gumbel(1) copula.
                indep_ent <- copmat == "indepcop"
                copmat[indep_ent] <- "gum"
                cparmat[indep_ent] <- list(1)
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
                res <- copreg::rVineTruncCondCDF(parvec = parvec,
                                         udat = dat,
                                         A = Aleaf,
                                         ntrunc = min(ntrunc, nrow(Aleaf) - 1),
                                         pcondmat = pcondmat,
                                         np = np)
                return(res)
            }
        }
    } else {
        ## Case 2: Vine cannot be subsetted to 'vbls'.
        ## Try removing as many variables as possible. Go through
        ##  the variables one at a time, each time a subset manages to happen.
        if (verbose)
            cat(paste0("Could not subset vine to variables ",
                       paste(vbls, collapse = ", "), ".\n"))
        want_rmvd <- setdiff(v, vbls)
        remaining <- want_rmvd
        while (length(remaining) > 0) {
            v_sub <- setdiff(v, remaining[1])
            rv_try <- subset(rv, v_sub)
            if (is.null(rv_try)) {
                remaining <- remaining[-1]
            } else {
                rv <- rv_try
                v <- v_sub
                want_rmvd <- setdiff(want_rmvd, remaining[1])
                remaining <- want_rmvd
            }
        }
        if (verbose)
            cat(paste0("Could only manage to subset to variables ",
                       paste(v, collapse = ", "), ".\n"))
        if (length(want_rmvd) + 1 > maxint)
            stop(paste0("Computing 'pcondrvine' requires more than ", maxint,
                        " integrals. Raise integration tolerance through 'maxint' argument."))
        if (verbose)
            cat(paste0("Obtaining density via ", length(want_rmvd)+1,
                       "-dimensional integration.\n"))
        ## Get density
        fX <- function(datvec) {  # Takes a full row of data.
            extradens <- function(uextra) {  # Accepts those variables to be removed.
                newvec <- datvec
                newvec[want_rmvd] <- uextra
                drvine(newvec, rv)
            }
            integrate.mv(extradens, rep(0, length(want_rmvd)), rep(1, length(want_rmvd)))
        }
        ## Get conditional distribution
        res <- apply(dat, 1, function(row) {
            dens <- function(xcond){  # Accepts uniform variable 'var'.
                x <- row
                x[var] <- xcond
                fX(x)
            }
            integrate.mv(dens, 0, row[var]) / integrate.mv(dens, 0, 1)
        })
        return(res)
    }
}
