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
