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

#' Fit bicop
#'
#' @param QYhat function of the previous
#' @param cpar list of copula parmeters you wish to fit to each family in \code{families}.
#' \code{NULL} if you don't want to specify any; otherwise, copulas for which
#' you don't want to specify should have a \code{NULL} entry; otherwise,
#' put a vector of at least partially specified parameters, with \code{NA}
#' going in place of parameters in that family that you don't know about.
#' @param taus Vector of quantile indices to score.
#' @import VineCopula
#' @export
fitbicop_cnqr <- function(u, v, uval, vval, QYhat,
                          families = c("indepcop", "bvncop","bvtcop","mtcj",
                                       "gum","frk","joe","bb1","bb7","bb8"),
                           cpar = NULL, taus = 1:10/11*0.1 + 0.9) {
    nfam <- length(families)
    ## Put cpar in list form if not already:
    if (is.null(cpar)) cpar <- rep(list(NULL), nfam)
    ## Initiate output
    valscore <- numeric(0)
    fittedcpars <- list()
    fittedcops <- character()  # We need this because the copula could be flipped.
    ## Fit each model
    for (i in 1:nfam) {
        cop <- families[i]
        ## Do we need to estimate any parameters? No, if the parameters are already
        ##  fully specified, or if there are no parameters.
        pre_cpar <- cpar[[i]]
        do_estimate <- TRUE
        if (!is.null(pre_cpar)) if (all(!is.na(pre_cpar))) {
            do_estimate <- FALSE
        }
        if (cop == "indepcop") {
            do_estimate <- FALSE
            pre_cpar <- numeric(0)
        }
        ## Were the parameters only partially specified?
        ## (CHECK)
        ## Estimate if need be.
        if (!do_estimate) {
            ## No need to estimate. Just get forecasts on the validation set.
            qcond <- get(paste0("qcond", cop))
            vhat_val <- sapply(taus, function(tau_) qcond(tau_, uval, pre_cpar))
        } else {
            ## We'll have to estimate. VineCopila::BiCopSelect has smart starting values,
            ##  so I'll use those to get starting values.
            startfit <- BiCopSelect(u, v, familyset = copname2num(cop)[[1]])
            ## Get proper copula rotation
            cop <- copnum2name(startfit$family)
            qcond <- get(paste0("qcond", cop))
            ## Extract starting values
            cparstart <- c(startfit$par, startfit$par2)
            if (cparstart[2] == 0) cparstart <- cparstart[1]
            ## Get forecasting family:
            ## ***NEED TO PUT IN PRE-SPECIFIED PORTION OF PARAMETER.
            forecast_fam <- function(u, theta) {
                sapply(taus, function(tau_) qcond(tau_, u, theta))
            }
            ## Choose theta to minimize score on training set:
            vhat_train <- function(theta) forecast_fam(u, theta)
            obj <- function(theta) scoreq(v, vhat_train(theta), taus)
            cparhat <- rnlm(obj, cparstart, cparspace(cop))$estimate
            ##
            ## Get forecast on validation set:
            vhat_val <- forecast_fam(uval, cparhat)
        }
        ## Get score from validation set:
        thisscore <- scoreq(vval, vhat_val, taus)
        ## Update outputs.
        valscore[i] <- thisscore
        fittedcpars[[i]] <- cparhat
        fittedcops[i] <- cop
    }
    ## Output
    list(cops = fittedcops, cpars = fittedcpars, scorehat = valscore)

}
