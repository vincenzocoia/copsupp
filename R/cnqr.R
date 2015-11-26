#' Fit a Non-linear Model using CNQR
#'
#' Regress a response against some predictors, using CNQR (composite
#' nonlinear quantile regression).
#' A non-linear model is found
#' using vine copulas -- see details.
#'
#' @param y Response vector.
#' @param xdat Matrix of predictor data. Variables are columns, observations
#' are in rows (which correspond to \code{y}).
#' @param tau Vector of quantile indices to regress on.
#' @param xmargs List of vectorized functions of marginal cdf's of the
#' predictor variables. Or, a single such function if they're all the same.
#' Or, leave \code{NULL} if you want to use the empirical distribution.
#' @param FY Vectorized function of the response cdf. Leave \code{NULL} if
#' you want to use the empirical distribution.
#' @param QY Vectorized function of the response quantile function. Leave \code{NULL} if
#' you want to use the empirical distribution.
#' @param ntruncX When a joint distribution of the response is sought using a vine,
#' what truncation level should be used (integer from 1 to \code{ncol(xdat)-1}).
#' Leave \code{NULL} for maximum truncation.
#' @param show_nlm After the cnqr optimization, should the \code{\link{nlm}}
#' results be displayed? \code{TRUE} if so.
#' @param familyset When choosing bivariate copula models to build the
#' distributions, what families should be used? Should be a vector of integer
#' codes used in \code{VineCopula::RVineCopSelect}.
#' @note Since this function is intended to be a "quick" regression,
#' you might encounter an error on the optimization portion through
#' \code{\link{nlm}}, especially when the gumbel copula is involved.
#'
#' Sorry about the matrix that gets printed whenever the function is run.
#' I can't seem to keep \code{VineCopula::RVineCopSelect} quiet.
#' @return A function that accepts the following arguments:
#'
#' \itemize{
#'      \item \code{$x}: A matrix of new observations, as in \code{xdat}. Could
#'      be a vector if only one observation.
#'      \item \code{$taunew}: A vector of quantile indices to evaluate the
#'      forecast quantile function at.
#' }
#'
#' It returns a matrix of the conditional quantiles, where columns correspond
#' to the quantile indices in \code{$tau}. Could be a vector if there's only
#' one observation in \code{x}.
#' @details Here is how the function works.
#'
#' \enumerate{
#'      \item If needed, univariate marginals are found using \code{\link{marginal}}.
#'      \item A model for the joint distribution of the predictors is fit
#'      using \code{\link{fitX}}.
#'      \item A model for the Bayesian Network linking the response to the
#'      predictors is found using \code{\link{fitBN}}.
#'      \item A family of forecasts are made for the
#'      data using \code{\link{pcondseq.vine}} and \code{\link{qcondBN}}.
#'      \item The optimal forecast is selected using \code{\link{nlm}}
#'      on an objective function found by \code{\link{cnqrobj}}. Starting values
#'      are obtained from the bivariate fits in the function previously used,
#'      \code{\link{fitBN}}.
#'      \item The optimal forecaster is extended to any new data, again using
#'      the functions \code{\link{pcondseq.vine}} and \code{\link{qcondBN}}.
#' }
#' @examples
#' ## Get some simulated data:
#' library(CopulaModel)
#' set.seed(73646)
#'
#' p <- 5
#' ntrunc <- p-1
#' n <- 50
#'
#' A0 <- truncvarray(Dvinearray(p), ntrunc)
#' copmat0 <- makeuppertri("frk", ntrunc, p, "")
#' cparmat0 <- makeuppertri(3, ntrunc, p)
#'
#' dat <- fvinesim(50, A0, copmat0, cparmat0)
#' dat <- qexp(dat)
#' y <- dat[, 1]
#' xdat <- dat[, -1]
#'
#' ## Get forecaster:
#' tau <- space_taus(10, tau_c = 0)
#' Qhat <- cnqr(y, xdat, tau = tau, verbose = TRUE)
#' Qhat(xdat[1, ])
#' Qhat(head(xdat), tau = c(0.9, 0.95, 0.99))
#' @export
cnqr <- function(y, xdat, tau = space_taus(10),
                 xmargs = NULL, FY = NULL, QY = NULL, ntruncX = NULL, verbose = FALSE,
                 familyset = c(1:10, 13, 14, 16:20, 23, 24, 26:30, 33, 34, 36:40)) {
    if (is.vector(xdat)) xdat <- matrix(xdat, ncol = 1)
    p <- ncol(xdat)
    n <- length(y)
    if (n != nrow(xdat)) stop("Number of observations do not match.")
    ## Fit marginals if need be.
    if (is.null(xmargs)) {
        xmargs <- apply(xdat, 2, function(col) marginal(col)$cdf)
    }
    if (length(xmargs) == 1) xmargs <- rep(list(xmargs), p)
    if (is.null(FY) | is.null(QY)) {
        ymarg <- marginal(y)
        if (is.null(FY)) FY <- ymarg$cdf
        if (is.null(QY)) QY <- ymarg$qf
    }
    ## Get truncation if need be.
    if (is.null(ntruncX)) ntruncX <- p-1
    ## Fit a model for predictors:
    if (verbose) cat("Fitting a model to the predictors...\n")
    yu <- FY(y)
    xdatu <- xdat
    for (col in 1:p) xdatu[, col] <- xmargs[[col]](xdat[, col])
    fitX. <- fitX(xdatu, ntrunc = ntruncX, familyset = familyset)
    if (verbose){
        print(fitX.)
        cat("Done.\n")
    }
    ## Fit a model for the Bayesian Network:
    if (verbose) cat("Fitting a model for the Bayesian Network...\n")
    fitBN. <- fitBN(yu, xdatu, familyset=familyset)
    #### Get lengths of parameter vectors
    len <- sapply(fitBN.$cparstart, length)
    if (verbose){
        print(fitBN.)
        cat("Done.\n")
    }
    ## Get sequences of conditional distributions.
    if (verbose) cat("Evaluating the sequential conditional predictor cdfs...\n")
    Fcond <- pcondseq.vine(fitBN.$xord, xdatu, rvinefit=fitX.,
                           .print = verbose, .print.cdfmethod = verbose)
    if (verbose) {
        cat("Done.\n")
        cat("Optimizing forecaster...\n")
    }
    ## Make predictions on the data. To speed up computation, if each parameter
    ##  is of length 1, there's no need to construct a list of parameters.
    if (all(len == 1)) {
        yhat <- function(cparvec) qcondBN(tau, fitBN.$cops, cparvec, Fcond, QY=QY)
    } else {
        yhat <- function(cparvec){
            cpar <- list()
            parnum <- 0
            for (i in 1:length(len)) {
                np <- len[i]
                cpar[[i]] <- cparvec[parnum + seq_len(np)]
                parnum <- parnum + np
            }
            qcondBN(tau, fitBN.$cops, cpar, Fcond, QY=QY)
        }
    }
    obj <- cnqrobj(y, yhat, tau)
    startval <- c(fitBN.$cparstart, recursive = TRUE)
    res <- nlm(obj, startval)
    cparfit <- res$estimate
    if (verbose){
        print(res)
        cat("Done.\n")
    }
    ## Now make a forecaster given new data
    function(x, taunew=tau) {
        if (is.vector(x)) x <- matrix(x, nrow = 1)
        if (ncol(x) != p) stop(paste("x should have", p, "variables."))
        Fcond <- pcondseq.vine(fitBN.$xord, x, rvinefit=fitX., FX = xmargs)
        qcondBN(taunew, fitBN.$cops, cparfit, Fcond, QY=QY)
    }
}
