#' Evaluate a function at limit points
#'
#' Sometimes evaluating a function in R at limit points does not give
#' the proper limit. \code{eval_lims} manually evaluates the
#' desired function values at those limit points, and evaluates the function
#' normally otherwise.
#'
#' @param fun Vectorized univariate function that you want to evaluate. Need not
#' evaluate to anything at points \code{replx}.
#' @param arg Vector of values to evaluate the function \code{fun} at.
#' @param replx Vector of values for which to manually evaluate the function
#' \code{fun} at.
#' @param replf Vector of values for which you want \code{replx} to evaluate to.
#' @param ... Other arguments to pass to the function \code{fun}.
#' Expecting these arguments to be vectors or lists. Any
#' vectorization in these arguments is preserved (see details to see how).
#' @return Vector of evaluations of the function \code{fun} at \code{arg}.
#' @note Any \code{NA}'s that appear in \code{arg} will be returned as
#' \code{NA}.
#' @details
#' To preserve vectorization over the \code{...} arguments, the vectors/lists
#' input here that match the length of \code{arg} are subsetted to match those
#' arguments in \code{arg} that are not the special values \code{replx}.
#' @examples
#' ## The function we want to evaluate, and some arguments:
#' fun <- function(x, alpha=1) exp(-1/x)/x^alpha
#' arg <- c(0:5, NA)
#'
#' ## Should have fun(0)=0, but we get NaN:
#' fun(arg)
#'
#' ## Manually evaluate:
#' eval_lims(fun, arg, replx=0, replf=0)
#'
#' ## Try other alpha values:
#' fun(arg, alpha=5)
#' eval_lims(fun, arg, replx=0, replf=0, alpha=5)
#' fun(arg, alpha=1:7)
#' eval_lims(fun, arg, replx=0, replf=0, alpha=1:7)
#'
#' ## NaN and numeric(0) work as arguments too:
#' eval_lims(fun, numeric(0), replx=0, replf=0)
#' eval_lims(fun, NaN, replx=0, replf=0)
#' @export
eval_lims <- function(fun, arg, replx, replf, ...) {
    ## Which of arg correspond to values in replx?
    ind <- lapply(replx, function(replx_) {
        which(arg == replx_)
    })
    allind <- c(ind, recursive=TRUE)
    if (length(allind) == 0) return(fun(arg, ...))
    ## Evaluate fun at the non-special arg values, and put them in the
    ##  results vector in the appropriate place:
    n <- length(arg)
    clean_arg <- arg[-allind]
    ellip <- list(...)
    clean_ellip <- lapply(ellip, function(vec) {
        if (length(vec) == n) return(vec[-allind]) else return(vec)
    })
    clean_f <- do.call(fun, c(list(clean_arg), clean_ellip)) #fun(clean_arg)
    res <- rep(NA, length(arg))
    res[-allind] <- clean_f
    ## Manually evaluate the function at the other points:
    for (i in seq_len(length(ind))) {
        ind_ <- ind[[i]]
        if (length(ind_) > 0) res[ind_] <- replf[i]
    }
    return(res)
}
