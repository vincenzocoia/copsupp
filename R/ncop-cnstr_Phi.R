#' Upper Incomplete Gamma Function
#'
#' The upper incomplete gamma function, defined as the integral
#' integrate(x^{k-1} exp(-x), x=t..infinity).
#'
#' @param k Indicated exponent in the integrand. Positive.
#' @param t Non-negative values of the lower limit of the integral.
#' @note This function is vectorized over either argument, but not both at once.
#' @export
igamma <- function(k, t) {
    gamma(k) * (1 - pgamma(t, k))
}



cnstr_Psiinv <- function(w, k, mxiter=20,eps=1.e-6,bd=5){
    ## Work with non-NA, non-1, non-0 values.
    NAs <- is.na(w)
    ones <- (w == 1)  # T/F. Has NA's too.
    whichones <- which(ones)
    zeroes <- (w == 0) # T/F. Has NA's too.
    whichzeroes <- which(zeroes)
    clean_w <- na.omit(w[!(ones | zeroes)])
    ## Compute gamma(k-1)
    gkm1 <- gamma(k-1)
    ## Go ahead with the algorithm
    if (length(clean_w) > 0){
        v <- clean_w
        tt=2/v # maybe a good start
        iter=0
        diff=1
        ## Begin Newton-Raphson algorithm
        while(iter<mxiter & max(abs(diff))>eps){
            igam0 <- igamma(k, tt)
            igam1 <- igamma(k-1, tt)
            g <- (k-1) * gkm1 - igam0 + tt*igam1 - tt*v*gkm1
            gp <- igam1 - v * gkm1
            diff=g/gp
            tt=tt-diff
            while(max(abs(diff))>bd | any(tt<=0))
            { diff=diff/2; tt=tt+diff }
            iter=iter+1
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
        res[whichones] <- 0
        res[whichzeroes] <- Inf
    } else {
        res <- tt
    }
    res
}

C2g1inv <- function(tau, u, k) {
    cnstr_Psi(k, qgamma(1-tau, k-1)/u)
}




u <- runif(1000)
tau <- runif(1000)
x <- qnorm(u)
y <- function(k) qnorm(C2g1inv(tau, u, k))
plot(x, y(2000))
