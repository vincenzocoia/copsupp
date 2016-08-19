## For the IG copula, tried using a "helper function", whose inverse if I had
##  I could evaluate C2g1inv (i.e. qcondigcop).

#' This function has problems when both theta and k are large.
ig_helper_inv <- function(w, theta, k, mxiter=20,eps=1.e-6,bd=5){
    ## Work with non-1 and non-0 values of w.
    ones <- (w == 1)  # T/F. Has NA's too.
    whichones <- which(ones)
    zeroes <- (w == 0) # T/F. Has NA's too.
    whichzeroes <- which(zeroes)
    clean_w <- na.omit(w[!(ones | zeroes)])
    ## which values are NA:
    NAs_th <- is.na(theta)
    NAs_w <- is.na(w)
    NAs <- unique(NAs_th, NAs_w)
    ## Go ahead with the algorithm
    if (length(clean_w) > 0){
        v <- clean_w
        ## Empirically, it looks like this tt is a good start:
        tt <- (1-v)^(-1/(theta+1))
        iter <- 0
        diff <- 1
        ## Begin Newton-Raphson algorithm
        while(iter<mxiter & max(abs(diff))>eps){
            ## Helpful quantities
            thlogt <- theta * log(tt)
            ## Evaluate functions
            g <- tt * (1-v) - pgamma(thlogt, k-1, lower.tail=FALSE)
            gp <- 1 - v + dgamma(thlogt, k-1) * theta / tt
            diff <- g/gp
            tt <- pmax(1, tt-diff)
            while(max(abs(diff))>bd | any(tt<=0))
            { diff <- diff/2; tt <- tt+diff }
            iter <- iter+1
            #cat(iter,diff,tt,"\n")
        }
    } else {
        tt <- numeric(0)
    }

    ## Set up vector to be returned (start off with NA's)
    res <- rep(NA, length(w))
    special_indices <- c(which(NAs), whichones, whichzeroes)
    if (length(special_indices > 0)){
        res[-special_indices] <- tt
        res[whichones] <- Inf
        res[whichzeroes] <- 1
    } else {
        res <- tt
    }
    res
}

#' @param cpar c(theta, k)
qcondigcop <- function(tau, u, cpar) {
    theta <- cpar[1]
    k <- cpar[2]
    inv <- ig_helper_inv(tau, theta * (1-u), k)
    1 - cnstr_H(inv, theta, k)
}

## Tried to invert cnstr_H using Newton-Raphson. Has problems when theta is large
##  and k isn't. I later tried changing the domain of cnstr_H from (1,oo) to
##  (-oo, oo) by doing a double-log transform (so that I could "open up" the
##  function near 1), but this still didn't work (resulted in all sorts of
##  NaN's and 0s for gp).
cnstr_Hinv <- function(w, theta, k, mxiter=20,eps=1.e-6,bd=5){
    ## Work with non-1 and non-0 values of w.
    ones <- (w == 1)  # T/F. Has NA's too.
    whichones <- which(ones)
    zeroes <- (w == 0) # T/F. Has NA's too.
    whichzeroes <- which(zeroes)
    clean_w <- na.omit(w[!(ones | zeroes)])
    ## which values are NA:
    NAs_th <- is.na(theta)
    NAs_w <- is.na(w)
    NAs <- unique(NAs_th, NAs_w)
    ## Go ahead with the algorithm
    if (length(clean_w) > 0){
        v <- clean_w
        ## Starting point (half-way between 1 and 1/v, since 1/v is the upper bound).
        tt <- (1/v - 1)/2 + 1
        iter <- 0
        diff <- 1
        ## Begin Newton-Raphson algorithm
        while(iter<mxiter & max(abs(diff))>eps){
            ## Helpful quantities
            logt <- log(tt)
            arg <- 1/(theta * logt)
            ## Evaluate functions
            g <- cnstr_Psi(arg, k) - tt * v
            gp <- - arg/logt * cnstr_DPsi(arg, k) - v
            ## Check that there are no NaN's which would result from Inf*0 in gp
            ##  (could use is.nan, but I'd rather ensure that the reason gp is NaN
            ##   is due to t being 1)
            gpnan <- tt==1
            if (any(gpnan)) {
                ## Replace the NaN's with the limit of gp, which depends on k.
                if (k > 2) {
                    gp[gpnan] <- -v[gpnan]
                } else if (k == 2) {
                    gp[gpnan] <- -theta/2 - v[gpnan]
                } else {
                    gp[gpnan] <- rep(-Inf, sum(gpnan))
                }
            }
            diff <- g/gp
            tt <- pmax(1, tt-diff)
            # cat("-------\n", iter,"\n",diff, "\n",tt,"\n")
            while(max(abs(diff))>bd | any(tt<=0))
            { diff <- diff/2; tt <- tt+diff }
            iter <- iter+1

        }
    } else {
        tt <- numeric(0)
    }

    ## Set up vector to be returned (start off with NA's)
    res <- rep(NA, length(w))
    special_indices <- c(which(NAs), whichones, whichzeroes)
    if (length(special_indices > 0)){
        res[-special_indices] <- tt
        res[whichones] <- 1
        res[whichzeroes] <- Inf
    } else {
        res <- tt
    }
    res
}




#' This one uses the inverse of the helper function.
#' @rdname igcop
#' @export
qcondigcop <- function(tau, u, cpar) {
    ## Get "parameters" useful for this function.
    theta <- cpar[1]
    k <- cpar[2]
    param <- (1 - u) * theta
    n_u <- length(u)
    n_tau <- length(tau)
    ## Was only one tau input? If so, repeat it to match the length of u.
    if (n_tau == 1) tau <- rep(tau, n_u)
    ## Which taus are 0 and 1? Treat them specially. u is allowed to be 1,
    ##  because igcop_helper_inv can handle theta=0. u=0 results in evaluating
    ##  its theta parameter at (1-u)*theta=theta, so no problem there.
    zeroes <- tau == 0  # Keeps NAs
    ones <- tau == 1  # Keeps NAs
    whichzeroes <- which(zeroes)  # NAs not included.
    whichones <- which(ones)  # NAs not included.
    ## Which are NA? This goes for both tau and u!
    NAs <- is.na(tau) | is.na(u)
    which_NA <- which(NAs)
    ## Only consider non-NA values of tau in (0,1). Even if this subsets
    ##  to an empty vector, the following procedure should still work.
    rmv <- union(which_NA, c(whichzeroes, whichones))
    if (length(rmv) > 0) {
        clean_tau <- tau[-rmv]  # Note tau[-integer(0)] = numeric(0)
    } else {
        clean_tau <- tau
    }
    ## If u is of length 1, then we can use vectorized property of
    ##  igcop_helper_inv().
    if (n_u == 1) {
        arg <- exp(exp(igcop_helper_inv(clean_tau, param, k)))
        res1 <- 1 - cnstr_H(arg, theta, k)
    } else {
        ## u (and param) is not of length 1. Suppose it's of length tau.
        ## Now subset param to match those of clean_tau.
        if (length(rmv) > 0) {
            clean_param <- param[-rmv]  # Note tau[-integer(0)] = numeric(0)
        } else {
            clean_param <- param
        }
        ## Get inverse helper function for each param value.
        arg <- mapply(function(tau_, param_) {
            exp(exp(igcop_helper_inv(tau_, param_, k)))
        }, clean_tau, clean_param)
        ## And the "clean" quantile values:
        res1 <- 1 - cnstr_H(arg, theta, k)
    }
    ## Lastly, put NAs where they belong, and replace tau=0 and tau=1
    ##  with quantiles of 0 and 1.
    if (length(rmv) == 0) {
        return(res1)
    } else {
        res <- rep(NA, length(tau))
        res[whichzeroes] <- 0
        res[whichones] <- 1
        res[-rmv] <- res1
        return(res)
    }
}

#' Get a grid for funding qcond.
#'
#' @param u Single numeric in [0,1]
#' @param theta Single numeric >0.
#' @note Since an interpolation is used to invert the function, and
#' there's a different function for each \code{theta}, this function
#' only allows one theta.
get_grid_qcond <- function(clean_tau, u, cpar, ngrid.min=1000) {
    ## Get grid on domain of phi:
    v <- (2*(1:ngrid.min) - 1) / (2*ngrid.min)
    fn <- pcondigcop(v, u, cpar)


    minw <- min(clean_tau)
    maxw <- max(clean_tau)
    wgridmin <- max(minw - 1/ngrid.min, minw/2)
    wgridmax <- min(maxw + 1/ngrid.min, (1+maxw)/2)
    wgrid <- seq(wgridmin, wgridmax, length.out=ngrid.min)
    tgrid <- log(k/theta*log(1/(1-wgrid)))
    ## Evaluate phi at that grid:
    fn <- igcop_helper(tgrid, theta, k)
    ## Does this cover the w values? If not, it ought to: expand the tgrid
    ##  so that it does. First, extend the tgrid to the left.
    fnlist <- list(fn)
    tlist <- list(tgrid)
    fnmin <- min(fn)
    fnmax <- max(fn)
    wpair <- fn[1:2]
    tpair <- tgrid[1:2]
    incr <- 1/ngrid.min
    i <- 1
    while(fnmin > minw) {  # Can do this because minw > 0 since clean_tau > 0.
        i <- i + 1
        ## Estimate most recent slope
        slope <- max(diff(wpair) / diff(tpair), 0.1)
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
        wnew <- seq(fnmax, wgridmax, length.out=ngrid.min)
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


#' This one inverts the pcond function directly.
#' @rdname igcop
#' @export
qcondigcop <- function(tau, u, cpar) {
    ## Get "parameters" useful for this function.
    theta <- cpar[1]
    k <- cpar[2]
    param <- (1 - u) * theta
    n_u <- length(u)
    n_tau <- length(tau)
    ## Was only one tau input? If so, repeat it to match the length of u.
    if (n_tau == 1) tau <- rep(tau, n_u)
    ## Which taus are 0 and 1? Treat them specially. u is allowed to be 1,
    ##  because igcop_helper_inv can handle theta=0. u=0 results in evaluating
    ##  its theta parameter at (1-u)*theta=theta, so no problem there.
    zeroes <- tau == 0  # Keeps NAs
    ones <- tau == 1  # Keeps NAs
    whichzeroes <- which(zeroes)  # NAs not included.
    whichones <- which(ones)  # NAs not included.
    ## Which are NA? This goes for both tau and u!
    NAs <- is.na(tau) | is.na(u)
    which_NA <- which(NAs)
    ## Only consider non-NA values of tau in (0,1). Even if this subsets
    ##  to an empty vector, the following procedure should still work.
    rmv <- union(which_NA, c(whichzeroes, whichones))
    if (length(rmv) > 0) {
        clean_tau <- tau[-rmv]  # Note tau[-integer(0)] = numeric(0)
    } else {
        clean_tau <- tau
    }




    ## If u is of length 1, then we can use vectorized property of
    ##  igcop_helper_inv().
    if (n_u == 1) {
        arg <- exp(exp(igcop_helper_inv(clean_tau, param, k)))
        res1 <- 1 - cnstr_H(arg, theta, k)
    } else {
        ## u (and param) is not of length 1. Suppose it's of length tau.
        ## Now subset param to match those of clean_tau.
        if (length(rmv) > 0) {
            clean_param <- param[-rmv]  # Note tau[-integer(0)] = numeric(0)
        } else {
            clean_param <- param
        }
        ## Get inverse helper function for each param value.
        arg <- mapply(function(tau_, param_) {
            exp(exp(igcop_helper_inv(tau_, param_, k)))
        }, clean_tau, clean_param)
        ## And the "clean" quantile values:
        res1 <- 1 - cnstr_H(arg, theta, k)
    }
    ## Lastly, put NAs where they belong, and replace tau=0 and tau=1
    ##  with quantiles of 0 and 1.
    if (length(rmv) == 0) {
        return(res1)
    } else {
        res <- rep(NA, length(tau))
        res[whichzeroes] <- 0
        res[whichones] <- 1
        res[-rmv] <- res1
        return(res)
    }
}



qcondigcop <- function(tau, u, cpar, mxiter=20,eps=1.e-6,bd=5){
    ## Make tau and u of the same length
    n_tau <- length(tau)
    n_u <- length(u)
    n <- max(n_tau, n_u)
    if (n_tau == 1) tau <- rep(tau, n)
    if (n_u == 1) u <- rep(u, n)
    # ## Work with non-NA, non-1, non-0 values.
    # NAs_tau <- is.na(tau)
    # NAs_u <- is.na(u)
    # NAs <-
    # ones <- (tau == 1)  # T/F. Has NA's too.
    # whichones <- which(ones)
    # zeroes <- (tau == 0) # T/F. Has NA's too.
    # whichzeroes <- which(zeroes)
    # clean_tau <- na.omit(tau[!(ones | zeroes)])
    clean_tau <- tau
    ## Go ahead with the algorithm
    theta <- cpar[1]
    k <- cpar[2]
    if (length(clean_tau) > 0){
        ## Use halfway between independence and comonotonicity.
        tt <- (clean_tau + u)/2
        iter <- 0
        diff <- 1
        ## Begin Newton-Raphson algorithm
        while(iter<mxiter & max(abs(diff))>eps){
            ## Helpful quantities
            Hkinv <- cnstr_Hinv(1-tt, theta, k)
            arg <- theta * (1-u) * log(Hkinv)
            der <- cnstr_D1H(Hkinv, theta, k)
            dgam <- dgamma(arg, k-1) * theta * (1-u)
            ## Evaluate functions
            g <- pgamma(arg, k-1, lower.tail=FALSE) - (1-clean_tau) * Hkinv
            gp <- 1/der * (dgam/Hkinv + 1 - clean_tau)
            diff <- g/gp
            tt <- tt-diff
            tt <- pmax(pmin(tt, 1), 0)
            while(max(abs(diff))>bd | any(tt<=0))
            { diff <- diff/2; tt <- tt+diff }
            iter <- iter+1
            #cat(iter,diff,tt,"\n")
        }
    } else {
        tt <- numeric(0)
    }

    ## Set up vector to be returned (start off with NA's)
    res <- rep(NA, max(length(tau), length(u)))
    special_indices <- c(which(NAs), whichones, whichzeroes)
    if (length(special_indices > 0)){
        res[-special_indices] <- tt
        res[whichones] <- 1
        res[whichzeroes] <- 0
    } else {
        res <- tt
    }
    res
    return(tt)


}


#' Newton-Raphson, forget about special values of tau, transform scale to R.
qcondigcop <- function(tau, u, cpar, mxiter=20,eps=1.e-6,bd=5){
    ## Make tau and u of the same length
    n_tau <- length(tau)
    n_u <- length(u)
    n <- max(n_tau, n_u)
    if (n_tau == 1) tau <- rep(tau, n)
    if (n_u == 1) u <- rep(u, n)
    ## Go ahead with the algorithm
    theta <- cpar[1]
    k <- cpar[2]
    if (length(tau) > 0){
        ## Use halfway between independence and comonotonicity.
        vv <- (tau + u)/2
        tt <- qnorm(vv)
        iter <- 0
        diff <- 1
        ## Begin Newton-Raphson algorithm
        while(iter<mxiter & max(abs(diff))>eps){
            ## Helpful quantities
            Hkinv <- cnstr_Hinv(1-vv, theta, k)
            arg <- theta * (1-u) * log(Hkinv)
            der <- cnstr_D1H(Hkinv, theta, k)
            dgam <- dgamma(arg, k-1) * theta * (1-u)
            ## Evaluate functions
            g <- pgamma(arg, k-1, lower.tail=FALSE) - (1-tau) * Hkinv
            gp <- 1/der * (dgam/Hkinv + 1 - tau) * dnorm(tt)
            diff <- g/gp
            tt <- tt-diff
            while(max(abs(diff))>bd)  # | any(tt<=0))
            { diff <- diff/2; tt <- tt+diff }
            vv <- pnorm(tt)
            iter <- iter+1
            # cat(iter,diff,tt,"\n")
        }
    } else {
        tt <- numeric(0)
    }
    return(pnorm(tt))


}

#' Newton-Raphson, forget about special values of tau, regular scale.
qcondigcop <- function(tau, u, cpar, mxiter=20,eps=1.e-6,bd=5){
    ## Make tau and u of the same length
    n_tau <- length(tau)
    n_u <- length(u)
    n <- max(n_tau, n_u)
    if (n_tau == 1) tau <- rep(tau, n)
    if (n_u == 1) u <- rep(u, n)
    ## Go ahead with the algorithm
    theta <- cpar[1]
    k <- cpar[2]
    if (length(tau) > 0){
        ## Use halfway between independence and comonotonicity.
        tt <- (tau + u)/2
        iter <- 0
        diff <- 1
        ## Begin Newton-Raphson algorithm
        while(iter<mxiter & max(abs(diff))>eps){
            ## Helpful quantities
            Hkinv <- cnstr_Hinv(1-tt, theta, k)
            arg <- theta * (1-u) * log(Hkinv)
            der <- cnstr_D1H(Hkinv, theta, k)
            dgam <- dgamma(arg, k-1) * theta * (1-u)
            ## Evaluate functions
            g <- pgamma(arg, k-1, lower.tail=FALSE) - (1-tau) * Hkinv
            gp <- 1/der * (dgam/Hkinv + 1 - tau)
            diff <- g/gp
            tt <- tt-diff
            while(max(abs(diff))>bd | any(tt<=0) | any(tt>=1))
            { diff <- diff/2; tt <- tt+diff }
            iter <- iter+1
            # cat(iter,diff,tt,"\n")
            print(qplot(g, gp) + labs(title=iter))
        }
    } else {
        tt <- numeric(0)
    }
    return(tt)


}
