#' Fit a new layer
#'
#' Choose and fit copula models on a new layer. The edge ("array column") must be
#' pre-specified. Intended for internal use.
#'
#' @param dat Data matrix with Uniform margins.
#' @param basevine Object of type "rvine" of the already-fit base vine for which
#' the new layer is to be applied.
#' @param edges Vector -- new column of vine array (with node appearing first)
#' @param cops Vector or list of pre-specified
#' copula families for each edge. Put \code{NA}
#' to leave the edge unspecified. \code{NULL} for
#' fully unspecified. You're allowed to put more than one family
#' as candidates.
#' @param cpars Pre-specified copula parameters corresponding to some of the
#' specified copulas in \code{cops}. Put \code{NA} in place of parameters to
#' leave them unspecified. \code{NULL} for fully unspecified.
#' @param families Vector of candidate copula family names for those that are
#' \code{NA} or \code{NULL}.
#' @import VineCopula
#' @note Expecting smart input. So, ensure that \code{edges} has length at
#' least 2, and that edges[-1] are variables in \code{basevine}, and that
#' \code{cpars} are only specified when there's only one copula family to
#' choose from.
#' @details Edges are fit so that edges[1] is the "V" variable. So, copulas
#' are fit to (edges[2], edges[1]), then (edges[3],
#' edges[1]) | edges[2], etc. That's because when computing edges[1]|others,
#' "pcond" can be used instead of "pcond12".
#' @return List of fitted \code{$cops} and \code{$cpars}.
#' @export
fitlayer <- function(dat, basevine, edges, cops = NULL, cpars = NULL,
                     families = c("indepcop", "bvncop","bvtcop","mtcj","gum",
                                  "frk","joe","bb1","bb7","bb8")) {
    d <- length(edges)
    if (is.null(cops)) cops <- rep(NA, d-1)
    if (is.null(cpars)) cpars <- rep(list(NA), d-1)
    if (!is.list(cops)) cops <- as.list(cops)
    if (!is.list(cpars)) cpars <- as.list(cpars)
    startpars <- cpars
    ## Go through edges and choose the best copula models. Do so by adding
    ##  variables one at a time.
    # for (j in 2:d) for (i in 1:(j-1)) {
    u <- dat[, edges[2]]
    v <- dat[, edges[1]]
    for (i in 2:length(edges) - 1) {
        ## Get candidate copula families
        if (any(is.na(cops[[i]]))) cand <- families else cand <- cops[[i]]
        ## Fit edge if cpar is not fully specified.
        if (all(is.na(cpars[[i]]))) {
            thisfit <- fitbicop_lh(u, v, families = cand)
        } else {
            thisfit <- fitbicop_lh(u, v, families = cand, cpar = cpars[[i]])
        }
        ## This edge is now fit. Update and get new (u,v) variables if there's still more.
        cops[[i]] <- thisfit$cop
        cpars[[i]] <- thisfit$cpar
        if (i+1 < length(edges)) {
            pcond <- get(paste0("pcond", thisfit$cop))
            v <- pcond(v, u, cpars[[i]])
            u <- pcondrvine(dat, basevine, cond = edges[i+2], vbls = edges[2:(i+2)])
        }
    }
    ## Output fits
    list(cops = c(cops, recursive = TRUE), cpars = cpars)
}

#' Temporary Title
#'
#' @param cops List with vector entries of the copula families to try fitting
#' for each edge.
#' @param cpars NOT FUNCTIONAL YET.
#' List with entries being the list of parameters you want to enforce
#' on the copula models. \code{NULL} for no restrictions; \code{NULL} entry
#' for no restrictions on that edge; \code{NA} in a parameter vector
#' means that parameter needs fitting.
#' @param w NOT FUNCTIONAL YET.
#' @param QY Vectorized function of the fitted quantile function of Y.
#' @param families Vector of copula families.
#' For those edges in \code{cop} that don't have copula families
#' specially selected (i.e. have \code{NULL} entries), these families are used.
#' @examples
#' set.seed(123)
#' library(CopulaModel)
#' library(cnqr)
#' library(ggplot2)
#'
#' ## Get data
#' rv <- rvine(AtoG(Dvinearray(5)), "frk", 3)
#' dattr <- rrvine(200, rv)
#' datval <- rrvine(200, rv)
#'
#' ## Fit CNQR:
#' fit1 <- fitlayer_cnqr(list(tr=dattr, val=datval), edges = 5:1,
#'                       rv, QY = identity, stoponbest = FALSE)
#' fit2 <- fitlayer_cnqr(list(tr=dattr, val=datval), edges = 5:1,
#'                       rv, QY = identity, stoponbest = FALSE,
#'                       cops = list("frk", NULL, NULL, NULL))
#'
#' ## Compare the two fits
#' #### Copulas:
#' fit1$cops
#' fit2$cops
#' #### Score:
#' plotdat <- data.frame(fit = rep(c("fit1", "fit2"), each = 5),
#'                       npred = rep(0:4, 2),
#'                       score = c(fit1$scores, fit2$scores))
#' ggplot(plotdat, aes(npred, score, group = fit)) +
#'     geom_point() +
#'     geom_line(aes(colour = fit))
#' @export
fitlayer_cnqr <- function(dats, edges, basevine, QY, tauset = space_taus(10),
                          stoponbest = TRUE, showplots = TRUE,
                          w = function(u) 1, cops = NULL, cpars = NULL,
                          families = c("indepcop", "bvncop","bvtcop","mtcj","gum",
                                       "frk","joe","bb1")) {
    library(ggplot2)
    ## 1. Modify and get info from input
    dattr <- dats$tr
    datval <- dats$val
    ytr <- dattr[, edges[1]]
    yval <- datval[, edges[1]]
    ## Fill-in cops if not done already.
    nfam <- length(families)
    nedge <- length(edges) - 1  # Number of edges
    if (is.null(cops)) cops <- rep(list(NULL), nedge)
    cops <- as.list(lapply(cops, function(fams){
        if (is.null(fams)) families else fams
    }))
    ## Make cpars a list if it's not already.
#     if (!is.list(cpars)) {
#         cpars <- rep(list(NULL), nedge)
#     }
    ## 2. Get sequential predictor cdfs
    ucondtr <- pcondseq(dattr, edges[-1], basevine)
    ucondval <- pcondseq(datval, edges[-1], basevine)
    ## 3. Fit each edge, one at a time.
    #### Initiate response variables
    vcondtr <- dattr[, edges[1]]
    vcondval <- datval[, edges[1]]
    #### The previous quantile function (or "most recent" one) should be a
    ####   function of tau (a vector, or a matrix with different taus
    ####   corresponding to observations), conddat (a matrix of uniform
    ####   conditional sequential predictors), and theta (a vector of
    ####   parameters).
    QYmat <- function(tau, conddat, theta) {
        if (is.vector(tau)) {
            QY_eval <- QY(tau)
            return(matrix(QY_eval, nrow = nrow(conddat), ncol = length(tau), byrow = TRUE))
        } else {
            return(t(apply(tau, 1, function(taurow) QY(taurow))))
        }
    }
    QYgX <- QYmat
    theta_prev_start <- numeric(0)
    paramspace <- function(theta) TRUE
    #### Initiate final fits
    fit_cops <- character(0)
    fit_scores <- numeric(0)
    fit_cpars <- list()
    fit_ncpars <- integer(0)
    #### Score the marginal
    yhatval <- QYgX(tauset, ucondval, numeric(0))
    fit_scores[1] <- scoreq(yval, yhatval, tau = tauset)
    #### Now loop over edges.
    for (i in 1:nedge) {
        ## BEGIN i'th EDGE. GOAL: Find a copula model.
        print(paste0("Edge #", i))
        edge_cops <- cops[[i]]
        edge_scores <- numeric(0)
        edge_cpars <- list()
        edge_Thetas <- list()
        ## (A) Fit and score each copula:
        for (j in 1:length(edge_cops)) {
            thiscop <- edge_cops[j]
            print(thiscop)
            ## Get starting values and copula rotation
            initfit <- fitbicop_lh(ucondtr[, i], vcondtr,
                                   families = thiscop)
            thiscop <- initfit$cop
            edge_cops[j] <- thiscop  # Because the copula may have rotated.
            cparstart <- initfit$cpar
            ncpar <- length(cparstart)
            Thetastart <- c(cparstart, theta_prev_start)
            thisqcond <- get(paste0("qcond", thiscop))
            ## Get a function that gets predictions from parameters
            yhat_fam <- function(conddat, Theta) {
                adj_tauset <- sapply(tauset, function(tau_) {
                    thisqcond(tau_, conddat[, i], Theta[seq_len(ncpar)])
                })
                theta <- Theta[ncpar + seq_len(length(Theta) - ncpar)]
                QYgX(adj_tauset, conddat, theta)
            }
            ## Objective function
            obj <- function(Theta)
                scoreq(ytr, yhat_fam(ucondtr, Theta), tau = tauset)
            ## Get parameter estimate:
            if (ncpar == 0) {
                thisTheta <- numeric(0)
            } else {
                ## Full parameter space
                thisparamspace <- function(Theta) {
                    theta <- Theta[ncpar + seq_len(length(Theta) - ncpar)]
                    cparspace(thiscop)(Theta[seq_len(ncpar)]) &
                        paramspace(theta)
                }
                ## Minimize objective function
                thisTheta <- rnlm(obj, Thetastart, thisparamspace)$estimate
            }
            edge_Thetas[[j]] <- thisTheta
            edge_cpars[[j]] <- thisTheta[seq_len(ncpar)]
            ## Score on the validation set
            thisyhatval <- yhat_fam(ucondval, thisTheta)
            edge_scores[j] <- scoreq(yval, thisyhatval, tau = tauset)
        }
        ## (B) Print results of this edge, if asked.
        if (showplots) {
            ## Scores:
            print(plotcop_score(edge_cops, edge_scores, fit_scores[i]) +
                      labs(title = paste0("Validation Scores for Edge #", i)))
            ## Normal scores plots:
            print(plotcop_qcurve(ucondval[, i], vcondval,
                                 cops = edge_cops[order(edge_scores)],
                                 cpars = edge_cpars[order(edge_scores)],
                                 tauset = tauset) +
                      labs(title = paste0("Normal Scores Validation Plots of",
                                          " Candidate Models: Edge #", i),
                           x = paste(edges[i+1], "given",
                                     paste(edges[1+seq_len(i-1)], collapse=", ")),
                           y = paste0("Response (variable ", edges[1],
                                      ") given ",
                                      paste(edges[1+seq_len(i-1)], collapse=", "))))
        }
        ## (C) Choose the best copula.
        edge_best <- which(edge_scores == min(edge_scores))[1] # ([1] needed b.c. RVineCopSelect.)
        new_score <- edge_scores[edge_best]
        new_cop <- edge_cops[edge_best]
        new_cpar <- edge_cpars[[edge_best]]
        new_Theta <- edge_Thetas[[edge_best]]
        ## Is this score better than before? If not, stop if asked.
        if (new_score > fit_scores[i] & stoponbest) {
            ## Adding this predictor does not improve the score.
            ##  Remove this variable and downstream variables from
            ##  the procedure, and break the loop.
            ## (i) Reduce sequential conditional predictors
            ucondtr <- ucondtr[, seq_len(i-1)]
            ucondval <- ucondval[, seq_len(i-1)]
            if (i-1 == 1) {
                ucondtr <- matrix(ucondtr, ncol = 1)
                ucondval <- matrix(ucondval, ncol = 1)
            }
            ## (ii) Reduce "edges"
            edges <- edges[1:i]
            break
        }
        ## Add the best fit to the layer:
        fit_cops[i] <- new_cop
        fit_cpars[[i]] <- new_cpar
        fit_ncpars[i] <- length(new_cpar)
        fit_scores[i+1] <- new_score
        ## How much better are we doing?
        if (showplots) {
            print(qplot(0:i, fit_scores, geom = c("point", "line"),
                        xlab = "Number of Predictors", ylab = "Score",
                        main = "CNQR score adjustment (Validation Set)"))
        }
        ## (D) Get ready for the next edge:
        if (i < nedge) {
            #### Condition the response again.
            pcond <- get(paste0("pcond", new_cop))
            pcondfit <- function(v, u) pcond(v, u, new_cpar)
            vcondtr <- pcondfit(vcondtr, ucondtr[, i])
            vcondval <- pcondfit(vcondval, ucondval[, i])
            #### Get updated parameter space
            theta_prev_start <- new_Theta
            ncpar_new <- length(new_cpar)
            i_still <- i
            paramspace <- function(Theta) {
                res <- TRUE
                for (e in i_still) {
                    ncpar <- fit_ncpars[e]
                    cpar <- Theta[seq_len(ncpar)]
                    Theta <- Theta[ncpar + seq_len(length(Theta) - ncpar)]
                    cop <- fit_cops[e]
                    res <- res & cparspace(cop)(cpar)
                }
                res
            }
            #### Get up-to-date quantile function
            QYgX <- function(tau, conddat, Theta) {
                ndat <- nrow(conddat)
                if (is.vector(tau))
                    tau <- matrix(tau, nrow = ndat, ncol = length(tau),
                                  byrow = TRUE)
                ## Move through the chain of qconds:
                for (e in i_still:1) {
                    cpar <- Theta[1:fit_ncpars[e]]
                    Theta <- Theta[-(1:fit_ncpars[e])]
                    qcond <- get(paste0("qcond", fit_cops[e]))
                    qcondfit <- function(p, u) qcond(p, u, fit_cpars[[e]])
                    tau <- apply(tau, 2, function(taucol){
                        qcondfit(taucol, conddat[, e])
                    })
                }
                ## Put on the marginal:
                QYmat(tau)
            }
        }
    }
    ## Output
    list(xord = edges[-1],
         cops = fit_cops,
         cpars = fit_cpars,
         basevine = basevine,
         QY = QY,
         useq = list(tr = ucondtr,
                     val = ucondval),
         scores = fit_scores)
}

#
#
#     ## Initiate output
#     valscore <- numeric(0)
#     fittedcpars <- list()
#     fittedcops <- character()  # We need this because the copula could be flipped.
#     ## Fit each model
#     for (i in 1:nfam) {
#         cop <- families[i]
#         ## Do we need to estimate any parameters? No, if the parameters are already
#         ##  fully specified, or if there are no parameters.
#         pre_cpar <- cpar[[i]]
#         do_estimate <- TRUE
#         if (!is.null(pre_cpar)) if (all(!is.na(pre_cpar))) {
#             do_estimate <- FALSE
#         }
#         if (cop == "indepcop") {
#             do_estimate <- FALSE
#             pre_cpar <- numeric(0)
#         }
#         ## Were the parameters only partially specified?
#         ## (CHECK)
#         ## Estimate if need be.
#         if (!do_estimate) {
#             ## No need to estimate. Just get forecasts on the validation set.
#             qcond <- get(paste0("qcond", cop))
#             vhat_val <- sapply(taus, function(tau_) qcond(tau_, uval, pre_cpar))
#         } else {
#             ## We'll have to estimate. VineCopila::BiCopSelect has smart starting values,
#             ##  so I'll use those to get starting values.
#             startfit <- BiCopSelect(u, v, familyset = copname2num(cop)[[1]])
#             ## Get proper copula rotation
#             cop <- copnum2name(startfit$family)
#             qcond <- get(paste0("qcond", cop))
#             ## Extract starting values
#             cparstart <- c(startfit$par, startfit$par2)
#             if (cparstart[2] == 0) cparstart <- cparstart[1]
#             ## Get forecasting family:
#             ## ***NEED TO PUT IN PRE-SPECIFIED PORTION OF PARAMETER.
#             forecast_fam <- function(u, theta) {
#                 sapply(taus, function(tau_) qcond(tau_, u, theta))
#             }
#             ## Choose theta to minimize score on training set:
#             vhat_train <- function(theta) forecast_fam(u, theta)
#             obj <- function(theta) scoreq(v, vhat_train(theta), taus)
#             cparhat <- rnlm(obj, cparstart, cparspace(cop))$estimate
#             ##
#             ## Get forecast on validation set:
#             vhat_val <- forecast_fam(uval, cparhat)
#         }
#         ## Get score from validation set:
#         thisscore <- scoreq(vval, vhat_val, taus)
#         ## Update outputs.
#         valscore[i] <- thisscore
#         fittedcpars[[i]] <- cparhat
#         fittedcops[i] <- cop
#     }
#     ## Output
#     list(cops = fittedcops, cpars = fittedcpars, scorehat = valscore)
# }
