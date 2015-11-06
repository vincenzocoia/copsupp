## Functions for the new copula construction
## C(u,v) = u + v - 1 + (1-u)*phi_{theta*(1-u)}(phiinv_{theta}(1-v))
#library(CopulaModel)

#' Construction ('phi') function for the new copula
#' 
#' phi[theta](t) = (1-t^(-theta)) / (t*theta*log(t))
#' 
#' @param theta Real number for the "shape" parameter of the function. Can't
#' be a vector.
#' @param t Vector of real numbers >= 1 to evaluate the function at.
phi <- function(theta, t){
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
  if (theta == 0){
    res[whichnotones] <- 1/t[whichnotones]
  } else {
    tt <- t[whichnotones]
    res[whichnotones] <- (1 - tt ^ (-theta)) / (tt * theta * log(tt)) 
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
phiinv <- function(theta, w, mxiter=20,eps=1.e-6,bd=5){
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
        ## "g" should be zero if tt is the solution.
        tth <- tt ^ (-theta)
        lt <- log(tt)
        g <- 1 - tth - theta * v * tt * lt
        ## Derivative of g wrt tt
        gp <- theta * (tth/tt - v*(lt+1))
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

#' Bivariate new copula -- cdf and density
#' 
#' \code{pnew} is the distribution function; \code{dnew} is the density;
#' \code{logdnew} is the log of the density.
#' 
#' @param u value in interval 0,1; could be a vector
#' @param v value in interval 0,1; could be a vector
#' @param cpar copula parameter (single number >= 0)
#' 
#' @details Either at least one of \code{u} or \code{v} should have length=1, or 
#' \code{u} and \code{v} should have the same length.
#' @rdname pnew
#' @export
pnew <- function(u, v, cpar){
  theprob <- function(one_u, one_v, one_cpar){
    one_u + one_v - 1 + (1-one_u) * 
      phi(one_cpar*(1-one_u), phiinv(one_cpar, 1-one_v))
  }
  res <- NA
  if (length(v) == 1){
    res <- sapply(u, function(uu) theprob(uu, v, cpar))
  } else {
    if (length(u) == 1){
      res <- sapply(v, function(vv) theprob(u, vv, cpar))
    }
    if (length(u) == length(v)){
      res <- sapply(1:length(u), function(i) theprob(u[i], v[i], cpar))
    }
  }
  res
}

#' @rdname pnew
#' @export
logdnew <- function(u, v, cpar){
  if (length(u) == 1 & length(v) > 1){
    u <- rep(u, length(v))
  }
  if (length(u) > 1 & length(v) == 1){
    v <- rep(v, length(u))
  }
  if (length(v) != length(u)){
    warning("u and v are not inputted with proper lengths. Returning NA.")
    return(NA)
  }
  alpha <- 1+cpar*(1-u)
  phiin <- phiinv(cpar, 1-v)
  log(alpha) - (alpha+1)*log(phiin) - log(-phip(cpar, phiin))
}

#' @rdname pnew
#' @export
dnew <- function(u, v, cpar){
  exp(logdnew(u, v, cpar))
}

#' Conditional distributions of bivariate new copula (\code{\link{pnew}}).
#' 
#' \code{pcondnew} and \code{pcondnew21} are the same functions -- they're
#' the cdf of V|U. \code{pcondnew12} is the cdf of U|V. \code{qcondnew} is
#' the quantile function of V|U ('U|V' not supported at this point, since
#' we only really care about V|U).
#' 
#' @param u value in interval 0,1; could be a vector
#' @param v value in interval 0,1; could be a vector
#' @param tau quantile indices in [0,1]; could be a vector.
#' @param cpar copula parameter (single number >= 0)
#' @details The first two arguments should have the same length. If not,
#' one of them should have length 1. 
#' @note \code{pcondnew21} is only included for completeness, since 
#' \code{pcondnew12} is introduced as well (being a non-symmetric copula).
#' @rdname pcondnew
#' @export
pcondnew <- function(v, u, cpar){
  theprob <- function(one_v, one_u, one_cpar){
    1 - phiinv(one_cpar, (1-one_v)) ^ (-one_cpar*(1-one_u)-1)
  }
  res <- NA
  if (length(v) == 1){
    res <- sapply(u, function(uu) theprob(v, uu, cpar))
  } else {
    if (length(u) == 1){
      res <- sapply(v, function(vv) theprob(vv, u, cpar))
    }
    if (length(u) == length(v)){
      res <- sapply(1:length(u), function(i) theprob(v[i], u[i], cpar))
    }
  }
  res
}

#' @rdname pcondnew
#' @export
pcondnew21 <- function(v, u, cpar)
  pcondnew(v, u, cpar)

#' @rdname pcondnew
#' @export
qcondnew <- function(tau, u, cpar){
  if (length(tau) == 1 & length(u) > 1){
    tau <- rep(tau, length(u))
  }
  if (length(tau) > 1 & length(u) == 1){
    u <- rep(u, length(tau))
  }
  if (length(tau) != length(u)){
    return(NA)
  }
  if (cpar == 0){
    return(tau)
  } else {
    alpha <- cpar*(1-u)+1
    return(1 - phi(cpar, (1-tau) ^ (-1/alpha)))
    ## Using this instead fixed a bug once:
    #a <- 1 - u + 1/cpar
    #t <- 1/(1-tau)
    #xi <- 1 / (cpar*(1-u) + 1)
    #return(1 - a * (1 - t^(-1/a)) / (t^xi * log(t)))
  }
}

#' @rdname pcondnew
#' @export
pcondnew12 <- function(u, v, cpar){
  theprob <- function(one_u, one_v, one_cpar){
    ub <- 1 - one_u
    vb <- 1 - one_v
    phiinvb <- phiinv(one_cpar, vb)
    1 - ub * phip(one_cpar * ub, phiinvb) / phip(one_cpar, phiinvb)
  }
  res <- NA
  if (length(v) == 1){
    res <- sapply(u, function(uu) theprob(uu, v, cpar))
  } else {
    if (length(u) == 1){
      res <- sapply(v, function(vv) theprob(u, vv, cpar))
    }
    if (length(u) == length(v)){
      res <- sapply(1:length(u), function(i) theprob(u[i], v[i], cpar))
    }
  }
  res
}

#' Simulate from the new copula family
#' 
#' @param n Positive integer: the number of observations to generate
#' @param cpar The value of the copula parameter (non-negative)
#' @return Returns an \code{n}x2 matrix. Each row is an independent draw
#' from the specified new copula.
#' @note This function is not vectorized.
#' @examples 
#' rnew(n = 10, cpar = 4.5)
#' @export
rnew <- function(n, cpar) {
  u <- runif(n)
  tau <- runif(n)
  v <- qcondnew(tau, u, cpar)
  matrix(c(u, v), ncol = 2, byrow = FALSE)
}

#' Kendall's tau of new copula family
#' 
#' Converts copula parameter of new copula to kendall's tau
#' 
#' @param cpar Vector of copula paremeters (>=0)
#' @return Vector of kendall tau values for the corresponding parameter values.
#' @details A wrapper for the function \code{\link{ktau}} in CopulaModel
#' package, using details of the new copula.
#' @import CopulaModel
#' @export
ktaunew <- function(cpar)
  sapply(cpar, function(theta) 
    ktau(theta, pcond12 = pcondnew12, pcond21 = pcondnew21, 
         dcop = dnew, pcop = pnew))
