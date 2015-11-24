#' Make CNQR objective function
#'
#' Makes the CNQR objective function for you.
#'
#' @param y A vector of observed responses/outcomes.
#' @param yhat_fun A function of the parameter vector that returns a matrix
#' with \code{length(y)} rows of forecasts corresponding to \code{y}.
#' Columns should be the quantile forecasts corresponding to \code{taus}.
#' @param taus A vector specifying the quantile indices that are being
#' predicted.
#' @param g A vectorized non-decreasing function for which to transform
#' both response/outcome and forecast in the estimation.
#' @return Returns a function that takes parameter values and returns
#' the value of the objective function.
#' @note You don't really need to use this function, since it's almost
#' identical to the \code{\link{scoreq}} function.
#' @seealso \code{\link{space_taus}}
#' @export
cnqrobj <- function(y, yhat_fun, taus, g = identity) {
    function(theta) scoreq(y, yhat_fun(theta), tau = taus, g = g)
}

#' Get equally spaced quantile indices.
#'
#' Obtain \code{K} equally spaced quantile indices in the open interval
#' \code{(tau_c, 1)}.
#'
#' @param K Integer. Number of quantile indices to obtain.
#' @param tau_c Value in the open interval (0,1) to take as the lower
#' quantile index bound.
#' @return A vector of length \code{K} with the quantile indices.
#' @export
space_taus <- function(K, tau_c = 0.9) {
    1:K/(K+1) * (1 - tau_c) + tau_c
}
