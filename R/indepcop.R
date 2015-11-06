#' Conditional Independence Copula
#' 
#' These are just the identity functions with respect to the first argument.
#' @rdname indepcopcond
#' @export
pcondindepcop <- function(v,u,cpar=0){
  v
}

#' @rdname indepcopcond
#' @export
qcondindepcop <- function(tau,u,cpar=0){
  tau
}

#' @rdname indepcopcond
#' @export
qcondindep <- function(tau,u,cpar=0){
  tau
}