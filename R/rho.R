#' Asymmetric absolute deviation function
#' 
#' The loss function that defines the tau^th-quantile;
#' rho(tau, u) = (tau - I(u<0))*u. Vectorised over \code{u}.
#' 
#' @param tau Number in [0,1] representing the probability/order
#' of the quantile.
#' @param u Value to evaluate the loss function at. Vectorised.
#' @export
rho <- function(tau, u){
  (tau - (u<0)) * u
}