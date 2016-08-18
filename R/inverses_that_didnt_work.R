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
