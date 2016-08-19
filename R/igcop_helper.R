igcop_helper <- function(t, theta, k) {
    1 - pgamma(theta*exp(t), k-1, lower.tail=FALSE) * exp(-exp(t))
}

#' Inverse of the helper function
#'
#' @param theta Single numeric >0.
#' @note Since an interpolation is used to invert the function, and
#' there's a different function for each \code{theta}, this function
#' only allows one theta.
get_grid_helper <- function(clean_w, theta, k, ngrid.min=1000) {
    ## Get grid on domain of phi:
    minw <- min(clean_w)
    maxw <- max(clean_w)
    wgridmin <- max(minw - 1/ngrid.min, minw/2)
    wgridmax <- min(maxw + 1/ngrid.min, (1+maxw)/2)
    wgrid <- seq(wgridmin, wgridmax, length.out=ngrid.min)
    tgrid <- log(k/theta*log(1/(1-wgrid)))
    ## Evaluate phi at that grid:
    fn <- igcop_helper(tgrid, theta, k)
    ## Does this cover the w values? If not, it ought to: expand the tgrid
    ##  so that it does.
    fnlist <- list(fn)
    tlist <- list(tgrid)
    fnmin <- min(fn)
    fnmax <- max(fn)
    wpair <- fn[1:2]
    tpair <- tgrid[1:2]
    incr <- 1/ngrid.min
    i <- 1
    while(fnmin > minw) {  # Can do this because minw > 0 (w=0 is treated specially)
        i <- i + 1
        ## Estimate most recent slope
        slope <- max(diff(wpair) / diff(tpair), 0.01)
        ## Fill-in the gap between the desired lower wgrid value and the actual
        ##  lower grid value (fnmin) with a new grid of values. But fnmin
        ##  can be 1 sometimes, so take wgridmax in that case.
        wnew <- seq(wgridmin, min(fnmin, wgridmax), length.out=ngrid.min)
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
    tgrid <- tgrid[i:1]
    fn <- fn[i:1]
    ## Continue if coverage above does not occur:
    fnlist <- list(fn)
    tlist <- list(tgrid)
    fnmin <- min(fn)
    fnmax <- max(fn)
    wpair <- tail(fn[[i]], 2)
    tpair <- tail(tgrid[[i]], 2)
    while(fnmax < maxw) {
        i <- i + 1
        ## Estimate most recent slope
        slope <- max(diff(wpair) / diff(tpair), 0.01)
        ## Fill-in the gap between the desired upper wgrid value and the actual
        ##  upper grid value (fnmax) with a new grid of values.
        wnew <- seq(fnmax, wgridmax, length.out=ngrid.min)
        ## Interpolate using the slope to find the t values that would evaluate
        ##  to wnew on the linear approximation.
        tnew <- tpair[2] + (fnmax - wnew) / slope
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

#' #' "Swallow" w.
#' #' @param n.init How many points should the algorithm begin with?
#' #' @param n.add How many points should the algorithm add to each iteration?
#' get_grid_helper <- function(clean_w, theta, k, n.add=1000, n.init=1000) {
#'     ## Get initial grid on the domain of 'helper' (note, clean_w might be length 1):
#'     approx_inv <- function(w) log(k/theta*log(1/(1-w)))
#'     wgrid_init <- ((1:n.init)*2-1)/(2*n.init)
#'     tgrid <- list(approx_inv(wgrid_init))
#'     ## Evaluate the function at the tgrid
#'     fn <- list(igcop_helper(tgrid[[1]], theta, k))
#'     ## Does the t grid cover the spread of inputted w's? We'll need these values
#'     ##  to find out:
#'     minw <- min(clean_w)
#'     maxw <- max(clean_w)
#'     minfn <- min(fn[[1]])
#'     maxfn <- max(fn[[1]])
#'     ## --Algorithm 1--
#'     ## Set up the algorithm to keep casting values to the *left* of tgrid until
#'     ##  the tgrid covers all w values.
#'     wpair <- fn[[1]][1:2]
#'     tpair <- tgrid[[1]][1:2]
#'     i <- 1
#'     while (minfn > minw) {
#'         i <- i + 1
#'         ## Estimate most recent slope
#'         slope <- max(diff(wpair) / diff(tpair), 0.01)
#'         ## Fill-in the gap between the desired lower wgrid value and the actual
#'         ##  lower grid value (minfn) with a new grid of values.
#'         wnew <- ((1:n.add)*2 - 2) / (2*n.add) * minfn # spread btwn (0, minfn).
#'         ## Interpolate using the slope to find the t values that would evaluate
#'         ##  to wnew on the linear approximation.
#'         tnew <- tpair[1] - (minfn - wnew) / slope
#'         ## Evaluate the actual igcop_helper function at these tnew values
#'         fnnew <- igcop_helper(tnew, theta, k)
#'         ## Calculate the minimum, and append these new values to the grid.
#'         wpair <- fnnew[1:2]
#'         tpair <- tnew[1:2]
#'         minfn <- wpair[1]
#'         fn[[i]] <- fnnew
#'         tgrid[[i]] <- tnew
#'     }
#'     ## Combine the fn and t values
#'     tgrid <- c(tgrid[i:1], recursive=TRUE)
#'     fn <- c(fn[i:1], recursive=TRUE)
#'     matrix(c(tgrid, fn), ncol=2)
#' }


#' @param theta Single numeric >0.
#' @param silent Logical; should the message output by \code{pcinterpolate()}
#' be silenced?
igcop_helper_inv <- function(w, theta, k, ngrid.min=1000, silent=TRUE) {
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
    grideval <- get_grid_helper(clean_w, theta, k)
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
