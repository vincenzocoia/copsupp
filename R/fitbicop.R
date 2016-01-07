#' Fit a Bivariate Copula Model
#'
#' These functions choose and fit the best copula model for data with
#' uniform margins. Intended for internal use.
#'
#' @param u Vector of Unif(0,1) data.
#' @param v Vector of Unif(0,1) data, same length as \code{u}.
#' @param families Vector of copula families to try fitting. Rotations
#' of each family are automatically considered.
#' @param method Either "aic" or "bic" (case doesn't matter)
#' @note This function is not that smart -- it assumes not all cpars are
#' specified, and if some are specified, that the copula model is known.
#'
#' I'm using BiCopSelect() just because it has smart starting points.
#' @return A list with entries \code{$cop} (the best-fitting copula family),
#' \code{$cpar} (the fitted copula parameters), and \code{$aic} and \code{$bic}.
#' @examples
#' library(CopulaModel)
#' set.seed(13)
#' u <- runif(100)
#' v <- qcondbvtcop(runif(100), u, c(0.9, 1))
#' fitbicop_lh(u, v)
#'
#' v <- qcondgum(runif(100), u, 4)
#' fitbicop_lh(u, 1-v)
#'
#' v <- runif(100)
#' fitbicop_lh(u, v)
#' @rdname fitbicop
#' @export
fitbicop_lh <- function(u, v, families = c("indepcop", "bvncop","bvtcop","mtcj",
                                             "gum","frk","joe","bb1","bb7","bb8"),
                        cpar = NULL, method = "aic") {
    ## Choose copula and get starting values
    #### Remove indepcop first. BiCopSelect doesn't like that.
    is_indepcop <- families == "indepcop"
    non_indep_cops <- families[!is_indepcop]
    startfit <- VineCopula::BiCopSelect(u, v, familyset = c(copname2num(non_indep_cops),
                                                            recursive = TRUE),
                                        selectioncrit = toupper(method))
    vinecopcpars <- c(startfit$par, startfit$par2)
    if (vinecopcpars[2] == 0) vinecopcpars <- vinecopcpars[-2]
    bestcop <- copnum2name(startfit$family)
    logdcop <- get(paste0("logd", bestcop))
    lastchar <- substring(bestcop, nchar(bestcop), nchar(bestcop))
    firstchars <- substring(bestcop, 1, nchar(bestcop)-1)
    if (exists(paste0("p", firstchars)) & lastchar %in% c("u", "v")) {
        vinecopcpars <- -vinecopcpars
    }
    ## If parts of the parameter are specified, re-do the fitting:
    if (!is.null(cpar)) {
        ## Which parameters are pre-specified?
        pre_spec <- !is.na(cpar)
        ## Get objective function
        obj <- function(theta) {
            fulltheta <- vinecopcpars # Dummy
            fulltheta[!pre_spec] <- theta
            fulltheta[pre_spec] <- cpar[pre_spec]
            -sum(logdcop(u, v, fulltheta))
        }
        finalfit <- nlm(obj, vinecopcpars[!pre_spec])
        cpar[!pre_spec] <- finalfit$estimate
    } else {
        cpar <- vinecopcpars
    }
    ## Get AIC or BIC. Note that AIC and BIC of indepcop is 0.
    nllh <- -sum(logdcop(u, v, cpar))
    npar <- length(cparspace(bestcop, F)[[1]])
    aic <- 2 * (npar + nllh)
    bic <- 2 * nllh + npar * log(length(u))
    if (tolower(method) == "aic") {
        criteria <- aic
    }
    if (tolower(method) == "bic") {
        criteria <- bic
    }
    ## If "indepcop" was selected as an option, check if it's AIC/BIC is better.
    if (any(is_indepcop)) if(criteria > 0) {
        ## Indepcop is better. Choose it.
        bestcop <- "indepcop"
        cpar <- integer(0)
        aic <- 0
        bic <- 0
    }
    ## Output results
    list(cop = bestcop, cpar = cpar, aic = aic, bic = bic)
}
