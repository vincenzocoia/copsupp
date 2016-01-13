#' Fit a Bivariate Copula Model
#'
#' These functions choose and fit the best copula model for data with
#' uniform margins. Intended for internal use. \code{fitbicop_qcond}
#' not done yet.
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
#' fitbicop_lh(u, v, families = "bvtcop", cpar = c(NA, NA))
#'
#' v <- qcondgum(runif(100), u, 4)
#' fitbicop_lh(u, 1-v)
#'
#' v <- runif(100)
#' fitbicop_lh(u, v)
#' fitbicop_lh(u, v, families = "indepcop")
#' @import VineCopula
#' @rdname fitbicop
#' @export
fitbicop_lh <- function(u, v, families = c("indepcop", "bvncop","bvtcop","mtcj",
                                             "gum","frk","joe","bb1","bb7","bb8"),
                        cpar = NULL, method = "aic") {
    ## Choose copula and get starting values
    #### Remove indepcop first. BiCopSelect doesn't like that.
    is_indepcop <- families == "indepcop"
    if (all(is_indepcop)) return(list(cop = "indepcop", cpar = integer(0), aic=0, bic=0))
    non_indep_cops <- families[!is_indepcop]
    #### Use BiCopSelect.
    startfit <- BiCopSelect(u, v, familyset = c(copname2num(non_indep_cops),
                                                            recursive = TRUE),
                                        selectioncrit = toupper(method))
    #### Extract chosen copula, and starting estimates.
    vinecopcpars <- c(startfit$par, startfit$par2)
    if (vinecopcpars[2] == 0) vinecopcpars <- vinecopcpars[-2]
    bestcop <- copnum2name(startfit$family)
    logdcop <- get(paste0("logd", bestcop))
    lastchar <- substring(bestcop, nchar(bestcop), nchar(bestcop))
    firstchars <- substring(bestcop, 1, nchar(bestcop)-1)
    if (exists(paste0("p", firstchars)) & lastchar %in% c("u", "v")) {
        vinecopcpars <- -vinecopcpars
    }
    ## Re-do the fitting, even if no parameters are specified. I'm not sure
    ##  that BiCopSelect chooses the MLE.
    if (is.null(cpar)) cpar <- rep(NA, length(cparspace(bestcop, F)[[1]]))
    #### Which parameters are pre-specified?
    pre_spec <- !is.na(cpar)
    #### Get and fit objective function
    obj <- function(theta) {
        fulltheta <- vinecopcpars # Dummy
        fulltheta[!pre_spec] <- theta
        fulltheta[pre_spec] <- cpar[pre_spec]
        if (!cparspace(bestcop)(fulltheta)) {
            99999999
        } else {
            -sum(logdcop(u, v, fulltheta))
        }
    }
    finalfit <- nlm(obj, vinecopcpars[!pre_spec])
    cpar[!pre_spec] <- finalfit$estimate
    ## Get AIC and BIC.
    nllh <- -sum(logdcop(u, v, cpar))
    npar <- length(cparspace(bestcop, F)[[1]]) - sum(pre_spec)
    aic <- 2 * (npar + nllh)
    bic <- 2 * nllh + npar * log(length(u))
    ## If "indepcop" was inputted as a candidate, check if it's AIC/BIC is
    ##  lower. Note that AIC and BIC of indepcop is 0.
    if (tolower(method) == "aic") {
        criteria <- aic
    }
    if (tolower(method) == "bic") {
        criteria <- bic
    }
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

#' @param taus Vector of quantile indices to score.
#' @rdname fitbicop
#' @export
fitbicop_qcond <- function(u, v, families = c("indepcop", "bvncop","bvtcop","mtcj",
                                              "gum","frk","joe","bb1","bb7","bb8"),
                           cpar = NULL, taus = 1:10/11*0.1 + 0.9) {

}
