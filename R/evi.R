#' Hill's estimator based on largest observations
#' 
#' From a univariate sample, returns the Hill estimate based on
#' the largest k+1 observations.
#' 
#' @param sample Numeric vector representing the univariate sample.
#' @param k Integer in [1,n] where \code{n=length(sample)}.
#' @param SEmethod Method used to estimate the standard error of the
#' estimate. If standard error is not required, put \code{NULL}. See details.
#' @return if \code{SEmethod=NULL}, returns a single numeric estimate of
#' the tail index xi. Otherwise, returns a 1x2 data frame of headers
#' \code{estimate} and \code{SE}.
#' @details The methods to estimate the standard error of the estimate is
#' currently \code{'norm'}, which uses the asymptotic Gaussian approximation
#' evaluated at the estimate of the tail index.
#' @seealso See the package 'smoothtail' for potentially more estimates.
#' @export
#' @examples
#' evi.Hill(rcauchy(1000), 100)
#' evi.Hill(rcauchy(1000), 100, NULL)
#' ## Can't compute if some values are non-positive:
#' evi.Hill(rnorm(1000) - 500, 100)
evi.Hill <- function(sample, k, SEmethod='norm'){
  n <- length(sample)
  ## Compute estimate
  ordered <- sort(sample)
  largest <- ordered[(n-k):n]
  if (largest[1] <= 0){
    warning("Cannot compute Hill Estimator; (n-k)th largest observation non-positive. Returning NA.")
    return(NA)
  }
  fraction <- largest/largest[1]
  xihat <- mean(log(fraction))
  ## Compute SE if required.
  if (!is.null(SEmethod)){
    ## Estimate SE
    if (SEmethod=='norm') {
      SE <- xihat/sqrt(k)
    }
    ## Output
    out <- data.frame(estimate = xihat, SE = SE)
  } else {
    out <- xihat
    names(out) <- "estimate"
  }
  out
}


#' Pickands Estimator based on largest observations
#' 
#' Computes the Pickands estimate of the extreme value index based
#' on a univariate sample, using the largest \code{4k+1} observations.
#' 
#' @param sample Numeric vector representing the univariate sample.
#' @param four_k Integer in [1,(n-1)] where \code{n=length(sample)}.
#' @param SEmethod Method used to estimate the standard error of the
#' estimate. If standard error is not required, put \code{NULL}. See details.
#' @return if \code{SEmethod=NULL}, returns a single numeric estimate of
#' the tail index xi. Otherwise, returns a 1x2 data frame of headers
#' \code{estimate} and \code{SE}.
#' @details The methods to estimate the standard error of the estimate is
#' currently \code{'norm'}, which uses the asymptotic Gaussian approximation
#' evaluated at the estimate of the tail index.
#' If \code{four_k} is not a multiple of 4, it will be reduced to
#' the nearest multiple of 4 that is less than \code{four_k}.
#' @seealso See the package 'smoothtail' for potentially more estimates.
#' @export
#' @examples
#' evi.Pickands(runif(1000), 50)
#' evi.Pickands(rcauchy(1000), 50)
#' evi.Pickands(runif(1000), 50, NULL)
evi.Pickands <- function(sample, four_k, SEmethod='norm'){
  ## First, reduce four_k to the next smallest multiple of 4, if it isn't
  ##  a multiple of 4 already
  k <- four_k - (four_k %% 4)
  ## Sample size
  n <- length(sample)
  ## Compute estimate
  ord <- sort(sample)
  xihat <- log((ord[n-k] - ord[n-2*k])/(ord[n-2*k] - ord[n-4*k])) / log(2)
  ## Compute SE if required.
  if (!is.null(SEmethod)){
    ## Estimate SE
    if (SEmethod=='norm') {
      SE <- xihat * sqrt(2^(2*xihat + 1) + 1) / log(4) / (2^xihat - 1) / sqrt(k)
    }
    ## Output
    out <- data.frame(estimate = xihat, SE = SE)
  } else {
    out <- xihat
    names(out) <- "estimate"
  }
  out
}