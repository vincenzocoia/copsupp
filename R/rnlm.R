#' Restricted \code{\link{nlm}}
#' 
#' A wrapper for \code{\link{nlm}} that puts a very large number wherever the
#' objective function does not exist.
#' 
#' @param f Function to minimize
#' @param p Starting point for the minimization
#' @param in.range A function that takes the
#'   argument of \code{f} and returns TRUE if that value of the 
#'   parameter is in the support of \code{f}, and returns FALSE if not.
#' @param lrg.number The large number that's put overtop of the objective
#' function outside of its support.
#' @param ign.error Sometimes \code{\link{nlm}} will throw an error (for example,
#' if it encounters numbers that are too big). Set this argument to \code{TRUE} to
#' have \code{NA} returned in the case of an error (and allow your code to keep
#' running afterwards). 
#' @param ... Other arguments to pass to \code{\link{nlm}}.
#' @export
rnlm <- function(f, p, in.range=NULL, lrg.number = 99999999, ign.error = FALSE, ...){
  ## Specify new objective function over the reals, based on the provided support.
  if (!is.null(in.range)) {
    ## If a restriction is specified:
    robj <- function(x){
      if (in.range(x)){
        return(f(x))
      } else {
        return(lrg.number)
      }
    }
  } else {
    ## Is a restriction is not specified:
    robj <- f
  }
  ## Minimize:
  #### First check if errors should be checked for or not.
  if (ign.error) {  # Yes -- we want to keep going if there's an error. 
    nlmres <- try(nlm(robj, p, ...))
    if (inherits(nlmres, "try-error")){
      ## There's been an error
      return(list(minimum = NA,
                  estimate = NA,
                  gradient = NA,
                  hessian = NA,
                  code = NA))
    } else {
      ## No error -- continue.
      return(nlmres)
    }
  } else {  # No -- let the error be thrown if there is one.
    return(nlm(robj, p, ...))
  }
}