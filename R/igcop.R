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

pcondigcop <- function(v, u, cpar) {
    theta <- cpar[1]
    k <- cpar[2]
    Hkinv <- cnstr_Hinv(1-v, theta, k)
    1 - pgamma(theta * (1-u) * log(Hkinv)) / Hkinv
}

qcondigcop <- function(tau, u, cpar) {
    ## Make tau and u of the same length
    n_tau <- length(tau)
    n_u <- length(u)
    n <- max(n_tau, n_u)
    if (n_tau == 1) tau <- rep(tau, n)
    if (n_u == 1) u <- rep(u, n)
    ## Use uniroot.
    sapply(seq_len(n), function(i) {
        f <- function(v) pcondigcop(v, u[i], cpar) - tau[i]
        # v0 <- 1 - cnstr_H((1-tau[i])^(-1/(theta+1)), theta, k)
        uniroot(f, c(0, 1))$root
    })
}
