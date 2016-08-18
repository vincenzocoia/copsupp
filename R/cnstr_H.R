cnstr_H <- function(t, theta, k) {
    cnstr_Psi(1/(theta * log(t)), k) / t
}

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





cnstr_Hinv <- function(w, theta, k, ngrid=1000) {
    ## Work with non-1 and non-0 values of w.
    ones <- (w == 1)  # T/F. Has NA's too.
    whichones <- which(ones)
    zeroes <- (w == 0) # T/F. Has NA's too.
    whichzeroes <- which(zeroes)
    clean_w <- na.omit(w[!(ones | zeroes)])
    ## which values are NA:
    NAs_th <- which(is.na(theta))
    NAs_w <- which(is.na(w))
    NAs <- unique(NAs_th, NAs_w)
    ## Main Idea: -----
    ## Get values to evaluate cnstr_H at, by using an approximation function
    ##  that can be inverted. We'll invert that function on a grid in (0,1)
    ##  to get a grid on the domain of H.
    ## It seems to be easier to find an approximation function by transforming
    ##  the domain from (1,oo) to (-oo,oo) by taking a log transform twice,
    ##  so we'll do that.
    ## ----------------
    ## Get grid on domain of H:
    wmin <- min(0.001, min(clean_w)/2)
    wmax <- max(0.999, (1+max(clean_w))/2)
    wgrid <- seq(wmin, wmax, length.out=ngrid)
    tgrid <- -log(-log(1-wgrid)) - log(1 + theta/k)
    ## Evaluate H at that grid:
    fn <- cnstr_H(exp(exp(tgrid)), theta, k)
    ## The evaluated cnstr_H become the domain values of cnstr_Hinv, and
    ##  the tgrid values become the evaluated cnstr_Hinv function.
    ## Use the pcinterpolate function to evaluate at w.
    der <- pcderiv(fn, tgrid)
    ## Evaluate at the requested NA-free w's:
    res1 <- pcinterpolate(fn, tgrid, der, clean_w)[, 1]
    ## Put back onto (1,oo) domain:
    res1 <- exp(exp(res1))
    ## Put in NA's where they were found:
    if (length(NAs) > 1) {
        res <- rep(NA, length(w))
        res[-NAs] <- res1
    } else {
        res <- res1
    }
    ## Manually evaluate the inverse function at zero and one:
    res[whichones] <- 1
    res[whichzeroes] <- Inf
    return(res)
}
