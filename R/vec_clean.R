#' 'Clean' a vector for function evaluation
#'
#'  Sometimes evaluating a function in R at limit points does not give
#'  what you want (such as a limit). This function manually inputs the
#'  desired function values at those points, and evaluates the function
#'  normally otherwise.
#'
#'  @param fun Vectorized univariate function that you want to evaluate. Need not
#'  evaluate to anything at points \code{replx}.
#'  @param arg Vector of values to evaluate the function \code{fun} at.
#'  @param replx Vector of values for which to manually evaluate the function
#'  \code{fun} at.
#'  @param replf Vector of values for which you want \code{replx} to evaluate to.
#'  @return Vector of evaluations of the function \code{fun} at \code{arg}.
#'  @note Any \code{NA}'s that appear in \code{arg} will be returned as
#'  \code{NA}.
#'  @examples
#'  ## The function we want to evaluate, and some arguments:
#'  fun <- function(x) exp(-1/x)/x
#'  arg <- c(0:5, NA)
#'
#'  ## Should have fun(0)=0, but we get NaN:
#'  fun(arg)
#'
#'  ## Manually evaluate:
#'  eval_clean(fun, arg, replx=0, replf=0)
#'
#'  ## NaN and numeric(0) work as arguments too:
#'  eval_clean(fun, numeric(0), replx=0, replf=0)
#'  eval_clean(fun, NaN, replx=0, replf=0)
#'  @export
eval_clean <- function(fun, arg, replx, replf) {
    ## Which of arg correspond to values in replx?
    ind <- lapply(replx, function(replx_) {
        which(arg == replx_)
    })
    allind <- c(ind, recursive=TRUE)
    if (length(allind) == 0) return(fun(arg))
    ## Evaluate fun at the non-special arg values, and put them in the
    ##  results vector in the appropriate place:
    clean_arg <- arg[-allind]
    clean_f <- fun(clean_arg)
    res <- rep(NA, length(arg))
    res[-allind] <- clean_f
    ## Manually evaluate the function at the other points:
    for (i in seq_len(length(ind))) {
        ind_ <- ind[[i]]
        if (length(ind_) > 0) res[ind_] <- replf[i]
    }
    return(res)
}
