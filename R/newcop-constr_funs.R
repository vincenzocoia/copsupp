#' Construction ('phi') function for the new copula
#'
#' phi[theta,k](t) = ell[k](theta * log(t)) / t
#'
#' @param theta Real number for the "shape" parameter of the function. Can't
#' be a vector.
#' @param t Vector of real numbers >= 1 to evaluate the function at.
#' @examples
#' cnstr(1:10, c(4,4))
#' cnstr(2:10, c(4,4))
#' cnstr(1, c(4,4))
#' cnstr(1:10, c(0,4))
#' cnstr(1:10, c(4,1))
#' ## As long as k doesn't go to 1 before theta goes to 0, we have:
#' cnstr(1:10, c(0,1))
#' @export
cnstr_H <- function(t, cpar){
    ## Which of the inputted t's =1?
    ones <- t == 1
    ## In order to accomodate possible NA's in the t vector, need to convert
    ##  the T/F vector "ones" to "which", which automatically does na.omit after.
    ## Indexes of t's that are 1:
    whichones <- which(ones)
    ## Index of t's that are !=1:
    whichnotones <- which(!ones)
    ## Set up
    res <- rep(NA, length(t))
    res[whichones] <- 1
    if (cpar[1] == 0){
        res[whichnotones] <- 1/t[whichnotones]
    } else {
        tt <- t[whichnotones]
        res[whichnotones] <- ell(cpar[1] * log(tt), cpar[2]) / tt
    }
    res
}

#' Inverse of the construction function -- old version with uniroot.
#'
#' Inverse of function 'phi'.
#' @param theta Real number for the "shape" parameter of the function. Can't
#' be a vector.
#' @param w Number in (0,1] to find the inverse at. Cannot be a vector.
phiinvold <- function(theta, w){
    ## Since phi is bounded above by 1/t, use that (shifted right by 0.5) as
    ##  the upper endpoint of the search interval.
    upper <- 1/w + 0.5
    uniroot(function(t) phi(theta, t) - w, c(1, upper))$root
}

#' Inverse of the construction function
#'
#' Uses 'newer' version suggested by Harry. It's possible that this
#' function could get stuck in an infinite loop, but it's never happened
#' before.
#'
#' @param theta Non-negative index for the construction function.
#' @param w Vector of values in in (0,1] to find the inverse at.
phiinv <- function(w, theta, k, mxiter=20,eps=1.e-6,bd=5){
    if (theta == 0) {
        res <- 1/w
    } else {
        ## Work with non-NA, non-1, non-0 values.
        NAs <- is.na(w)
        ones <- (w == 1)  # T/F. Has NA's too.
        whichones <- which(ones)
        zeroes <- (w == 0) # T/F. Has NA's too.
        whichzeroes <- which(zeroes)
        clean_w <- na.omit(w[!(ones | zeroes)])
        ## Go ahead with the algorithm
        if (length(clean_w) > 0){
            v <- clean_w
            tt=1/v # maybe a good start
            iter=0
            diff=1
            while(iter<mxiter & max(abs(diff))>eps){
                thlogt <- theta * log(tt)
                g <- ell(k, thlogt) - v * tt
                gp <- theta * ellp(k, thlogt) / tt - v
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
            res[whichones] <- 1
            res[whichzeroes] <- Inf
        } else {
            res <- tt
        }
    }
    res
}



#' Derivative of the construction function (with respect to argument)
#'
#' diff(phi[theta](t), t)
#'
#' @param theta Non-negative index for the phi function.
#' @param t Numeric values >=1 to evaluate the function at. Could be a vector.
phip <- function(theta, t){
    ones <- t == 1
    res <- rep(NA, length(t))
    res[ones] <- -(1+theta/2)
    tt <- t[!ones]
    res[!ones] <- (tt ^ (-theta-1) - (1+log(tt))*phi(theta, tt)) / (tt * log(tt))
    res
}
