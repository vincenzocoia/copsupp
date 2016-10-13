#' IG Copula helper function
#'
#' These functions are useful for computing the 2|1
#' distribution of the IG copula family. We have
#' C_{2|1}(v|u) = igcop_helper(log(log(x)), theta*(1-u), k), where
#' x=cnstr_Hinv(1-v, theta, k),
#' so that the quantile function can be written in terms of
#' \code{igcop_helper_inv}. The inverse is found using
#' interpolation, for which \code{\link{get_grid_helper}} finds a grid of
#' \code{t} values such that \code{igcop_helper} evaluates to points
#' that contain the \code{w}'s input into \code{igcop_helper_inv}.
#'
#' @param t Numeric vector to evaluate \code{igcop_helper} at. Real.
#' @param w Numeric vector to evaluate the inverse function at,
#' \code{igcop_helper_inv}. In [0,1].
#' @param theta Single numeric >0.
#' @param k Single numeric >1.
#' @param silent Logical; should the message output by \code{pcinterpolate()}
#' be silenced?
#' @note Since an interpolation is used to invert \code{igcop_helper},
#' which is a different function for each \code{theta},
#' \code{igcop_helper_inv} only allows one theta.
#' @return \code{igcop_helper} and \code{igcop_helper_inv} both return
#' vectors of length equal to the length of their first argument, representing
#' the function evaluations.
#' @rdname igcop_helper
#' @export
igcop_helper <- function(t, theta, k) {
    et <- exp(t)
    1 - pgamma(theta*et, k-1, lower.tail=FALSE) * exp(-et)
}

#' Evaluate the helper function at a grid
#'
#' Chooses a grid on the real line so that evaluating the function
#' \code{\link{igcop_helper}} at those points ends up spanning its range (0,1)
#' and encapsulating \code{w}. It's used so that
#' \code{\link{igcop_helper}} can be inverted by interpolation.
#'
#' @param w Vector of values in (0,1) for which the grid should cover. Should not
#' contain \code{NA}s.
#' @param theta Single numeric >0; parameter of \code{\link{igcop_helper}}
#' @param k Single numeric >1; parameter of \code{\link{igcop_helper}}
#' @param ngrid Positive integer; how many points at a time should the
#' interpolation grid be expanded by?
#' @details This function begins by obtaining a first guess at a grid whose
#' evaluated \code{\link{igcop_helper}} would encapsulate \code{w}. It
#' does this by inverting a rudimentary approximation function of
#' \code{\link{igcop_helper}},  1-exp(-theta/k*exp(t)) for t real. Then
#' if need be,
#' the grid is expanded by adding more points, chosen
#' by linear extrapolation based on the slopes
#' at either end of the grid.
#' @note This algorithm is somewhat rudimentary. It would be nice to upgrade it
#' so that it finds a more "focussed" grid that doesn't necessarily have to
#' span the entire range of the \code{\link{igcop_helper}} function, but just
#' the range of \code{w}. But if \code{w} is too narrow, the algorithm
#' is susceptible for being too sparse for \code{\link{igcop_helper}}
#' to be interpolated accurately. So \code{w} is currently expanded if it's not
#' "wide enough".
#' @return \code{get_grid_helper} returns a matrix with
#' 2 columns: column one contains values of \code{t}, and column two
#' contains \code{igcop_helper} evaluated at those \code{t} values (which
#' should "engulf" \code{clean_w}).
#' @examples
#' ## Grid for the theta=k=3 function:
#' foo <- get_grid_helper(0.4, 3, 3)
#' plot(foo[, 1], foo[, 2])
#'
#' foo <- get_grid_helper(1:99/100, 3, 3) # Same thing.
#' plot(foo[, 1], foo[, 2])
#'
#' ## Here's one where the grid needed expanding:
#' foo <- get_grid_helper(0.4, 0.1, 1.1)
#' nrow(foo)
#' plot(foo[, 1], foo[, 2])
#' @export
get_grid_helper <- function(w, theta, k, ngrid=1000) {
    ## Get starting grid on domain of igcop_helper:
    minw <- min(w)
    maxw <- max(w)
    wgridmin <- min(1/ngrid, minw/2)
    wgridmax <- max(1-1/ngrid, (1+maxw)/2)
    wgrid <- seq(wgridmin, wgridmax, length.out=ngrid)
    tgrid <- log(k/theta*log(1/(1-wgrid)))
    ## Evaluate igcop_helper at that grid:
    fn <- igcop_helper(tgrid, theta, k)
    ## Does this cover the w values? If not, it ought to: expand the tgrid
    ##  so that it does. First, extend the tgrid to the left.
    fnlist <- list(fn)
    tlist <- list(tgrid)
    fnmin <- min(fn)
    fnmax <- max(fn)
    wpair <- fn[1:2]
    tpair <- tgrid[1:2]
    incr <- 1/ngrid
    i <- 1
    while(fnmin > minw) {  # Can do this because minw > 0 since w > 0.
        i <- i + 1
        ## Estimate most recent slope
        slope <- max(diff(wpair) / diff(tpair), 0.1)
        ## Fill-in the gap between the desired lower wgrid value and the actual
        ##  lower grid value (fnmin) with a new grid of values. But fnmin
        ##  can be 1 sometimes, so take wgridmax in that case.
        wnew <- seq(wgridmin, min(fnmin, wgridmax), length.out=ngrid)
        ## Interpolate using the slope to find the t values that would evaluate
        ##  to wnew on the linear approximation.
        tnew <- tpair[1] - (fnmin - wnew) / slope
        ## Evaluate the actual igcop_helper function at these tnew values
        fnnew <- igcop_helper(tnew, theta, k)
        ## Calculate the minimum, and append these new values to the grid.
        wpair <- fnnew[1:2]
        tpair <- tnew[1:2]
        fnmin <- wpair[1]
        fnlist[[i]] <- fnnew
        tlist[[i]] <- tnew
    }
    ## Switch the order of t values and function values.
    tlist <- tlist[i:1]
    fnlist <- fnlist[i:1]
    ## Continue if we need tgrid to reach further right.
    wpair <- tail(fnlist[[i]], 2)
    tpair <- tail(tlist[[i]], 2)
    while(fnmax < maxw) {
        i <- i + 1
        ## Estimate most recent slope
        slope <- max(diff(wpair) / diff(tpair), 0.1)
        ## Fill-in the gap between the desired upper wgrid value and the actual
        ##  upper grid value (fnmax) with a new grid of values.
        wnew <- seq(fnmax, wgridmax, length.out=ngrid)
        ## Interpolate using the slope to find the t values that would evaluate
        ##  to wnew on the linear approximation.
        tnew <- tpair[2] + (wnew - fnmax) / slope
        ## Evaluate the actual igcop_helper function at these tnew values
        fnnew <- igcop_helper(tnew, theta, k)
        ## Calculate the minimum, and append these new values to the grid.
        wpair <- tail(fnnew, 2)
        tpair <- tail(tnew, 2)
        fnmax <- wpair[2]
        fnlist[[i]] <- fnnew
        tlist[[i]] <- tnew
    }
    tgrid <- c(tlist, recursive=TRUE)
    fn <- c(fnlist, recursive=TRUE)
    matrix(c(tgrid, fn), ncol=2)
}

#' @rdname igcop_helper
#' @import CopulaModel
#' @export
igcop_helper_inv <- function(w, theta, k, ngrid=1000, silent=TRUE) {
    if (length(w) == 0) return(numeric(0))
    if (theta == 0) return(log(log(1/(1-w))))
    ## Work with non-1 and non-0 values of w.
    ones <- (w == 1)  # T/F. Has NA's too.
    whichones <- which(ones)
    zeroes <- (w == 0) # T/F. Has NA's too.
    whichzeroes <- which(zeroes)
    clean_w <- na.omit(w[!(ones | zeroes)])
    ## which values are NA:
    NAs <- which(is.na(w))
    ## Main Idea: -----
    ## Get values to evaluate cnstr_H at, by using an approximation function
    ##  that can be inverted. We'll invert that function on a grid in (0,1)
    ##  to get a grid on the domain of H.
    ## It seems to be easier to find an approximation function by transforming
    ##  the domain from (1,oo) to (-oo,oo) by taking a log transform twice,
    ##  so we'll do that.
    ## ----------------
    ## Get grid and function evaluations:
    grideval <- get_grid_helper(clean_w, theta, k, ngrid = ngrid)
    tgrid <- grideval[, 1]
    fn <- grideval[, 2]
    ## The evaluated cnstr_H become the domain values of cnstr_Hinv, and
    ##  the tgrid values become the evaluated cnstr_Hinv function.
    ## Use the pcinterpolate function to evaluate at the requested
    ##  NA-free w's.
    if (silent) {
        ## The only way I can find to silent a cat() call is to use sink,
        ##  but send it to a non-existing file. Got answer from Stack
        ##  Overflow by cbielow:
        ## http://stackoverflow.com/questions/6177629/how-to-silence-the-output-from-this-r-package
        f <- file()
        sink(file=f)
        der <- pcderiv(fn, tgrid)
        res1 <- pcinterpolate(fn, tgrid, der, clean_w)[, 1]
        sink()
        close(f)
    } else {
        der <- pcderiv(fn, tgrid)
        res1 <- pcinterpolate(fn, tgrid, der, clean_w)[, 1]
    }

    ## Put in NA's where they were found:
    if (length(NAs) > 1) {
        res <- rep(NA, length(w))
        res[-NAs] <- res1
    } else {
        res <- res1
    }
    ## Manually evaluate the inverse function at zero and one:
    res[whichones] <- Inf
    res[whichzeroes] <- -Inf
    return(res)
}