#' Fit a Non-linear Model using CNQR
#'
#' Regress a response against some predictors. A non-linear model is found
#' using vine copulas -- see details.
cnqr <- function(y, xdat, xmargs = NULL, FY = NULL, QY = NULL, ntruncX = NULL,
                 familyset = c(1:10, 13, 14, 16:20, 23, 24, 26:30, 33, 34, 36:40)) {
    if (is.vector(xdat)) xdat <- matrix(xdat, ncol = 1)
    p <- ncol(xdat)
    n <- length(y)
    if (n != nrow(xdat)) stop("Number of observations do not match.")
    ## Fit marginals if need be.
    if (is.null(xmargs)) {
        xmargs <- apply(xdat, 2, function(col) marginal(col)$cdf)
    }
    if (is.null(FY) | is.null(QY)) {
        ymarg <- marginal(y)
        if (is.null(FY)) FY <- ymarg$cdf
        if (is.null(QY)) QY <- ymarg$qf
    }
    ## Get truncation if need be.
    if (is.null(ntruncX)) ntruncX <- p-1
    ## Fit a model for predictors:
    yu <- FY(y)
    xdatu <- xdat
    for (col in 1:p) xdatu[, col] <- xmargs[[col]](xdat[, col])
    fitX(xdat, ntrunc = ntruncX)

}
