#' Fit a Bayesian Network model
#'
#' Finds copula models for each pair of response and predictor in some order,
#' so that the pair is conditional on previous predictors.
#'
#' @param y Vector of response observations.
#' @param xdat Matrix of predictors. Columns are variables, rows are observations
#' (that correspond with entries of \code{y}).
#' @param ymarg Vectorized function of the cdf of the response.
#' @param xmargs List of vectorized functions of the cdfs of the univariate
#' predictors, in the order of the columns of \code{xdat}. Or, a single function
#' if they're all the same.
#' @param xord If you already know what order you want to pair the response
#' and the predictors, put that order here as a vector of "column numbers".
#' @param familyset A vector of integer codes of the copula families to try
#' fitting. See \code{VineCopula::BiCopSelect} for a full list.
#' @param ... Other arguments to pass to \code{VineCopula::BiCopSelect}.
#' @details This function first determines the order to pair up the response
#' and predictors in the order of highest partial correlation (by using
#' \code{\link{lm}}). Then the bivariate copula models are chosen and fitted
#' individually using \code{VineCopula::BiCopSelect}.
#'
#' By "pairing response and predictors in some order", I mean pairing (Y,X1),
#' (Y,X2)|X1, (Y,X3)|(X1,X2), etc, though not necessarily in that order.
#'
#' For the \code{familyset} argument, the default is almost all of
#' the families available. It just doesn't include the Tawn copula families.
#' @return A list with three entries:
#'
#' \itemize{
#'      \item \code{$xord}: The pairing order, as a vector of the column numbers
#'      of \code{xdat}.
#'      \item \code{$cops}: A vector of copula family names corresponding to
#'      the pairs in the order in \code{$xord}.
#'      \item \code{$cparstart}: A list of copula parameters corresponding to
#'      the families in \code{$cops}, which should be taken as starting values
#'      for an optimization (obtained by bivariate likelihood). If each entry
#'      of the list is a vector of length one, then a vector of those values
#'      is returned here instead.
#' }
#' @examples
#' ## Get some simulated data:
#' set.seed(6277)
#' p <- 5
#' A0 <- CopulaModel::Dvinearray(p)
#' copmat0 <- makeuppertri("frk", p-1, p, "")
#' cparmat0 <- makeuppertri(3, p-1, p)
#' dat <- fvinesim(100, A0, copmat0, cparmat0)
#' y <- dat[, 1]
#' xdat <- dat[, -1]
#'
#' ## Fit the model:
#' fit.BN(y, xdat)
#' @export
fit.BN <- function(y, xdat, ymarg = identity, xmargs = identity, xord = NULL,
                   familyset = c(1:10, 13, 14, 16:20, 23, 24, 26:30, 33, 34, 36:40), ...) {
    if (is.vector(xdat)) xdat <- matrix(xdat, ncol = 1)
    p <- ncol(xdat)
    if (length(xmargs) == 1) xmargs <- rep(list(xmargs), p)
    ## Uniformize the data
    for (col in 1:p) xdat[, col] <- xmargs[[col]](xdat[, col])
    y <- ymarg(y)
    ## Choose sequence of X's.
    if (is.null(xord)) {
        ynorm <- qnorm(y)
        xdatnorm <- qnorm(xdat)
        ## Just go along and see what has the highest (partial) correlations
        cor1 <- cor(ynorm, xdatnorm)
        xord <- which(cor1^2 == max(cor1^2))
        for (i in 1+seq_len(p-1)) { # 2:p if p>1.
            tocondn <- xdatnorm[, xord[i-1]]
            ynorm <- lm(ynorm ~ tocondn)$residuals
            xdatnorm <- apply(xdatnorm, 2, function(col) {
                if (!all(col==0)) lm(col ~ tocondn)$residuals else col
            })
            cors <- apply(xdatnorm, 2, function(col){
                if (!all(col==0)) cor(ynorm, col) else 0
            })
            xord <- c(xord, which(cors^2 == max(cors^2)))
            xdatnorm[, xord[i]] <- 0
        }
    }
    ## Now fit copula models and get starting points for parameter estimates
    cops <- character(0)
    cpars <- list()
    for (vbl in xord) {
        thisfit <- VineCopula::BiCopSelect(y, xdat[, vbl], familyset=familyset, ...)
        cops <- c(cops, copnum2name(thisfit$family))
        thiscpar <- numeric(0)
        if (thisfit$par != 0) thiscpar <- c(thiscpar, thisfit$par)
        if (thisfit$par2 != 0) thiscpar <- c(thiscpar, thisfit$par2)
        cpars <- c(cpars, list(thiscpar))
    }
    if (all(sapply(cpars, length) == 1)) cpars <- c(cpars, recursive = TRUE)
    ## Output:
    list(xord=xord, cops=cops, cparstart=cpars)
}
