#' Bivariate new copula -- cdf and density
#'
#' \code{pnew} is the distribution function; \code{dnew} is the density;
#' \code{logdnew} is the log of the density. This version is the 2-parameter
#' version, with cpar=c(theta, k).
#'
#' @param u value in interval 0,1; could be a vector
#' @param v value in interval 0,1; could be a vector
#' @param cpar copula parameters (theta>=0, k>=1).
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

#' Helper function for qcondnew
#'
#' This is the inverse of the function that maps t to
#' (1-pgamma(a*log(t),k-1))/t, which has domain [1, infinity) and range
#' (0,1].
#' @param w Value to evaluate the inverse at. Should be in [0,1] (with 0
#' returning \code{Inf}). Vectorized.
#' @param a,k Parameters, as above.
qcondnew_helper <- function(w, a, k, mxiter=20,eps=1.e-6,bd=5){
    if (a == 0) {
        return(1/w)
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
            tt=(1/v - 1)/2 + 1 # maybe a good start
            iter=0
            diff=1
            while(iter<mxiter & max(abs(diff))>eps){
                alogt <- a * log(tt)
                g <- 1 - pgamma(alogt, k-1) - v * tt
                gp <- - a * dgamma(alogt, k-1) / tt - v
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

#' #' @rdname pcondnew
#' #' @export
#' qcondnew <- function(tau, u, cpar){
#'     if (length(tau) == 1 & length(u) > 1){
#'         tau <- rep(tau, length(u))
#'     }
#'     if (length(tau) > 1 & length(u) == 1){
#'         u <- rep(u, length(tau))
#'     }
#'     if (length(tau) != length(u)){
#'         return(NA)
#'     }
#'     len <- length(tau)
#'     if (cpar[1] == 0){
#'         return(tau)
#'     } else {
#'         return(sapply(1:len, function(i)
#'             cnstr(qcondnew_helper(tau[i], cpar[1]*u[i], cpar[2]), cpar)))
#'     }
#' }



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
