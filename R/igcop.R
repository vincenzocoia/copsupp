







pcondigcop <- function(v, u, cpar) {
    theta <- cpar[1]
    k <- cpar[2]
    Hkinv <- cnstr_Hinv(1-v, theta, k)
    1 - pgamma(theta * (1-u) * log(Hkinv), k-1, lower.tail=FALSE) / Hkinv
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

#' @rdname cnstr_Psi
#' @export
qcondigcop <- function(tau, u, cpar, mxiter=20,eps=1.e-6,bd=5){
    ## Work with non-NA, non-1, non-0 values.
    NAs <- is.na(tau)
    ones <- (tau == 1)  # T/F. Has NA's too.
    whichones <- which(ones)
    zeroes <- (tau == 0) # T/F. Has NA's too.
    whichzeroes <- which(zeroes)
    clean_tau <- na.omit(tau[!(ones | zeroes)])
    ## Go ahead with the algorithm
    theta <- cpar[1]
    k <- cpar[2]
    if (length(clean_tau) > 0){
        ## Just use the identity function as a starting value.
        tt <- clean_tau
        iter <- 0
        diff <- 1
        ## Begin Newton-Raphson algorithm
        while(iter<mxiter & max(abs(diff))>eps){
            ## Helpful quantities
            Hkinv <- cnstr_Hinv(1-tt, theta, k)
            arg <- theta * (1-u) * log(Hkinv)
            der <- cnstr_DH(Hkinv, theta, k)
            dgam <- dgamma(arg, k-1) * theta * (1-u)
            ## Evaluate functions
            g <- dgam/der/Hkinv + (1-tau)/der
            gp <- 1/der * (dgam/Hkinv + 1 - tau)
            diff <- g/gp
            tt <- tt-diff
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
}
