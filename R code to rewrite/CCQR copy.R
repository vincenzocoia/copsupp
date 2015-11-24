## Make function to fit CNQR

#' Fit model via CNQR
#' 
#' Estimate parameters of a model specified for some upper range of
#' quantiles using CNQR estimation.
#' 
#' @param xdat A CanmoreR data object of covariates. Alternatively,
#' a covariate matrix (or vector if there's only one covariate).
#' @param ydat A CanmoreR data object of length 1 of response values.
#' Alternatively, a vector of response values (that match up with the rows
#' of the covariate matrix in \code{xdat}).
#' @param QYgX Quantile function of Y|X. Arguments: 
#' (single quantile index, x1 vector, x2 vector, vector of parameters). Should return
#' the quantile.
#' @param theta_init a vector of initial parameter value for optimization
#' @param tau_c Lower quantile index for which the \code{QYgX} model applies.
#' Not needed if \code{tau} is specified.
#' @param K Number of quantile curves to regress upon, with indices to be spread
#' evenly throughout \code{(tau_c, 1)}.
#' @param tau Vector of quantile indices to regress upon. If NULL, \code{K}
#' quantile indices will be chosen evenly spaced throughout \code{(tau_c, 1)}. 
#' Fit a model to bivariate data using CNQR.
#' 
#' @details Output
#' 
#' \code{$estimate}: Estimated model parameter(s)
#' 
#' \code{$objective}: Value of the objective function at the minimum
#' 
#' \code{$predict}: Prediction equation (tau, x) of the tau'th quantile curve
#' at x.
#' @export
cnqr <- function(xdat, ydat, QYgX, theta_init, tau_c=0.9, K=10, tau=NULL){
  if (length(ydat) > 1 & is.list(ydat)) {
    stop("More than one response was input into cnqr function; can only handle one.")
  }
  ## Get x data as a matrix, and y as a vector (but ensure the dates match up)
  if (is.list(ydat)){
    ## Data were entered as CanmoreR data objects
    dat <- MatchDates.slim2(c(ydat, xdat))
    y <- dat[, 1]
    X <- dat[, -1]
    n <- nrow(dat)
  } else {
    ## Data were supposedly entered as vectors and/or matrices
    X <- as.matrix(xdat)
    y <- ydat
    n <- length(y)
  }
  
  ## Determine quantile indices to regress upon
  if(is.null(tau)){
    ## Evenly spread K quantile indices throughout (tau_c, 1)
    tau <- (1:K)/(K+1)*(1-tau_c) + tau_c
  } else {
    if (any(tau>=1 | tau<=0)) stop('Inputted quantile indices not all in (0,1)')
    ## If vector of indices (tau) is specified, we just need to know its length.
    K <- length(tau)
  }
  
  ## Set up objective function
  ## Technique: sum the n x K matrix of rho_{tau[k]}(y[i] - QYgX(tau[k], x[i], theta))
  obj <- function(theta) {
    sum(sapply(tau, function(tau_) {
      ## Vector of quantiles:
      #QYgXvec <- apply(X, 1, function(row){
      #  QYgX(tau_, row[1], row[2], theta)
      #})
      QYgXvec <- QYgX(tau_, X[, 1], X[, 2], theta)
      ## 
      sum(rho(tau_, y - QYgXvec)) 
    }))
  }

  ## Minimize the function
  res <- nlm(obj, theta_init)
  names(res)[1] <- "objective"
  res
}

#' Fit bivariate (i.e. X,Y) model via CNQR
#' 
#' Estimate parameters of a bivariate model specified for some upper range of
#' quantiles using CNQR estimation.
#' 
#' @param x Vector of covariate values
#' @param y Vector of response values
#' @param QYgX Quantile function of Y|X. Arguments: 
#' (single quantile index, vector of X values to compute
#' \code{QYgX} for, vector of parameters)
#' @param theta_range If the model parameter is a real number and you want
#' to use \code{optimize} instead of \code{nlm}, specify a vector of length two
#' indicating the lower and upper bounds to search for solution.
#' @param theta_init If you have a multi-dimensional parameter, then
#' \code{nlm} needs to be used, so specify a vector of initial parameter value
#' (or specify this if you'd rather use \code{nlm})
#' @param tau_c Lower quantile index for which the \code{QYgX} model applies.
#' Not needed if \code{tau} is specified.
#' @param K Number of quantile curves to regress upon, with indices to be spread
#' evenly throughout \code{(tau_c, 1)}.
#' @param tau Vector of quantile indices to regress upon. If NULL, \code{K}
#' quantile indices will be chosen evenly spaced throughout \code{(tau_c, 1)}. 
#' Fit a model to bivariate data using CNQR.
#' @param plot_window If you want to see a plot of the objective function (so 
#' long as the parameter is 1D), specify the "window width" of the plot here.
#' Otherwise, keep it NULL.
#' 
#' @details Output
#' 
#' \code{$estimate}: Estimated model parameter(s)
#' 
#' \code{$objective}: Value of the objective function at the minimum
#' 
#' \code{$predict}: Prediction equation (tau, x) of the tau'th quantile curve
#' at x.
#' @import ggplot2
#' @export
cnqr_biv <- function(x, y, QYgX, theta_range, theta_init=NULL, tau_c=0.9, K=10, tau=NULL, plot_window=NULL){
  ## Look at the data
  n <- length(x)
  if(n != length(y)) 
    stop('x and y data are of different lengths.')
  which_na <- is.na(x) | is.na(y)
  if(any(which_na) > 0){
    x <- x[-which_na]
    y <- y[-which_na]
    warning(paste(sum(which_na), 
                   'observations with missing values removed from the data.'))
  }
  
  ## Determine quantile indices to regress upon
  if(is.null(tau)){
    ## Evenly spread K quantile indices throughout (tau_c, 1)
    tau <- (1:K)/(K+1)*(1-tau_c) + tau_c
  } else {
    if (any(tau>=1 | tau<=0)) stop('Inputted quantile indices not all in (0,1)')
    ## If vector of indices (tau) is specified, we just need to know its length.
    K <- length(tau)
  }
  
  ## Set up objective function
  #### n x K matrix of rho_{tau[k]}(y[i] - QYgX(tau[k], x[i], theta))
  v <- function(theta) sapply(tau, function(tau_) 
    rho(tau_, y - QYgX(tau_, x, theta)))
  #### Objective function
  obj <- function(theta) sapply(theta, function(theta_) sum(v(theta_)))
  
  ## Minimize the function
  if(is.null(theta_init)){
    res <- optimize(obj, theta_range)
    names(res)[1] <- "estimate"
  } else {
    res <- nlm(obj, theta_init)
    names(res)[1] <- "objective"
  }
  
  ## Plot the objective?
  if (!is.null(plot_window) & length(res$estimate) == 1) {
    estim <- res$estimate
    print(ggplot(data.frame(x = c(estim - plot_window/2, estim + plot_window/2)),
           aes(x)) + stat_function(fun = obj))
    #print(curve(obj, estim - plot_window/2, estim + plot_window/2))
  }
 
  res
}

#' Fit model via CCQR
#' 
#' Estimate parameters of a bivariate copula with specified margins, 
#' for some upper range of
#' quantiles using CNQR estimation.
#' For the estimation, response data are transformed to have Exp(1) margins.
#' 
#' @param x Vector of covariate values
#' @param y Vector of response values
#' @param qcond Quantile function of conditional copula, as in \code{qcond}
#' in the \code{CopulaModel} package.
#' @param theta_range If the model parameter is a real number and you want
#' to use \code{optimize} instead of \code{nlm}, specify a vector of length two
#' indicating the lower and upper bounds to search for solution.
#' @param theta_init If you have a multi-dimensional parameter, then
#' \code{nlm} needs to be used, so specify a vector of initial parameter value
#' (or specify this if you'd rather use \code{nlm})
#' @param FX Fully specified marginal distribution of the covariate,
#' vectorized. May cause problems if
#' it evaluates to 1 at some \code{x} data points.
#' Default if \code{NULL} is the empirical cdf, where the value at
#' each 'step' is between the step heights. 
#' @param FY Fully specified marginal distribution of the response,
#' vectorized. Will cause problems if
#' it evaluates to 1 at some \code{y} data points.
#' Default if \code{NULL} is the empirical cdf, where the value at
#' each 'step' is between the step heights. 
#' @param FYinv Fully specified quantile function of the response, vectorized.
#' If \code{NULL}, uses R's \code{\link{quantile}} function from the y data.
#' @param tau_c Lower quantile index for which the \code{QYgX} model applies.
#' Not needed if \code{tau} is specified.
#' @param K Number of quantile curves to regress upon, with indices to be spread
#' evenly throughout \code{(tau_c, 1)}.
#' @param tau Vector of quantile indices to regress upon. If NULL, \code{K}
#' quantile indices will be chosen evenly spaced throughout \code{(tau_c, 1)}. 
#' Fit a model to bivariate data using CNQR.
#' 
#' @details Output:
#' 
#' \code{$estimate}: Estimated model parameter(s)
#' 
#' \code{$objective}: Value of the objective function at the minimum
#' 
#' \code{$predict}: Prediction equation (tau, x) of the tau'th quantile curve
#' at x.
#' 
#' If \code{FY = NULL}, 
#' @export
ccqr <- function(x, y, qcond, theta_range, theta_init=NULL, FX = NULL,
                 FY = NULL, FYinv=NULL, tau_c=0.9, K=10, tau=NULL){
  n <- length(x)
  if(n != length(y)) 
    stop('x and y data are of different lengths.')
  ## Get edf marginals if NULL
  ## Note: Make a separate function that's 0.5/n lower than the edf
  ##        to apply to the data themselves. Otherwise it doesn't matter.
  if(is.null(FX)){
    FX <- function(xdata){
      val <- ecdf(x)(xdata)
      val[val == 1] <- 0.99999999999
      val[val == 0] <- 0.00000000001
      val
    }
    FXdata <- function(xdata) ecdf(x)(xdata) - 0.5/n
  } else {
    FXdata <- FX
  }
  if(is.null(FY)){
    FY <- function(ydata) ecdf(y)(ydata)
    FYdata <- function(ydata) FY(ydata) - 0.5/n
  } else {
    FYdata <- FY
  }
  if(is.null(FYinv)){
    FYinv <- function(tau) quantile(y, tau)
  }
  
  ## Transform data:
  u <- FXdata(x)
  yp <- qexp(FYdata(y))
  
  ## Get model:
  QYgX <- function(tau, u, cpar) qexp(qcond(tau, u, cpar))
  
  ## Fit:
  fit <- cnqr(u, yp, QYgX, 
              theta_range = theta_range, 
              theta_init = theta_init,
              tau_c = tau_c,
              K = K,
              tau = tau)
  estim <- fit$estimate
  
  ## Tack on the prediction equation too.
  fit$predict <- function(tau, x) FYinv(qcond(tau, FX(x), estim))
  
  fit
}