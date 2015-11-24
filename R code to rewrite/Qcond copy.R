#' Make Conditional Quantile Surface for 2-covariate simple vine
#' 
#' Constructs a function that is the quantile function 
#' of Y|U1,U2 (if \code{FX = NULL}) or
#' Y|X1,X2, where X1 and X2 are the covariates, and U1 and U2 are their
#' probability transformations. Vine structure: X1--X2--Y.
#' 
#' @param QY A vectorized function of the quantile function of the response.
#' @param cops copula names as in \code{\link{pcop}} in CopulaModel, like "frk".
#' Only need to specify one if all three copulas are the same.
#' See details for the order. 
#' @param FX Either \code{NULL} if the ouputted function is to take the
#' probability transforms \code{u1} and \code{u2} as arguments, or a list whose
#' two entries are the vectorized distribution functions of X1 and X2,
#' respectively.
#' @param numpars Vector of length three specifying the dimension of the
#' parameters that define the copula families (in the respective order).
#' @details
#' Specify copulas in this order: First for (X1,X2); Second for (X2,Y);
#' Third for (X1,Y)|X2. If the first is asymmetric, \code{Make_Qcond_2sv} 
#' assumes the cdf of 1|2 is named \code{pcondcop12}, 
#' where \code{cop} is the name of the copula, but if such an object does not exist,
#' then \code{pcondcop} is used instead.
#' @return
#' Returns a function that accepts the following arguments:
#' 
#' \code{tau_} A single quantile index
#' 
#' \code{u1} or \code{x1}: vector of the probability transformations of 
#' covariate X1 if \code{FX = NULL}, otherwise X1 itself.
#' 
#' \code{u2} or \code{x2}: vector of the probability transformations of 
#' covariate X2 if \code{FX = NULL}, otherwise X2 itself.
#' (Note that \code{u1/x1} and \code{u2/x2} must be of the same length.)
#' 
#' \code{cpars} vector of copula parameters corresponding to \code{cops}.
#' For copula families with parameter dimension !=1, just put the parameters
#' one after the other. (Need to do it this way to put into a numerical
#' optimizer).
#' @export
Make_Qcond_2sv <- function(QY, cops, FX = NULL, numpars = c(1,1,1)){
  if(length(cops) == 1) cops <- rep(cops, 3)
  ## Copula functions
  #### The "1|2" distribution of C1:
  pcond1_12full <- paste0("pcond", cops[1], "12")
  if (exists(pcond1_12full)) {
    pcond1_12 <- get(pcond1_12full)
  } else {
    pcond1_12 <- get(paste0("pcond", cops[1]))
  }
  #### The "2|1" quantile function of C2:
  qcond2_21 <- get(paste0("qcond", cops[2]))
  #### The "2|1" quantile function of C3:
  qcond3_21 <- get(paste0("qcond", cops[3]))
  ## The copula parameter indices from 'numpar':
  ind <- 1:sum(numpars)
  ind1 <- ind[seq(length.out = numpars[1])]
  ind <- setdiff(ind, ind1)
  ind2 <- ind[seq(length.out = numpars[2])]
  ind3 <- setdiff(ind, ind2)
  if (is.null(FX)){
    function(tau_, u1, u2, cpars){
      ## Parameters
      cpar1 <- cpars[ind1]
      cpar2 <- cpars[ind2]
      cpar3 <- cpars[ind3]
      ## Compute quantiles (in stages for easier debugging)
      comp1 <- pcond1_12(u1, u2, cpar1)
      comp2 <- qcond3_21(tau_, comp1, cpar3)
      comp3 <- qcond2_21(comp2, u2, cpar2)
      QY(comp3)
    }
  } else {
    function(tau_, x1, x2, cpars){
      ## Parameters
      cpar1 <- cpars[[1]]
      cpar2 <- cpars[[2]]
      cpar3 <- cpars[[3]]
      ## Probability transforms
      u1 <- FX[[1]](x1)
      u2 <- FX[[2]](x2)
      ## Compute quantiles (in stages for easier debugging)
      comp1 <- pcond1_12(u1, u2, cpar1)
      comp2 <- qcond3_21(tau_, comp1, cpar3)
      comp3 <- qcond2_21(comp2, u2, cpar2)
      QY(comp3)
    }
  }
}


#' Make CNQR objective function
#'
#' Function that builds an objective function. So, the output is
#' the objective function.
#'
#' @param Qcond The conditional quantile surface. Should have arguments
#' as in the output of \code{\link{Make_Qcond_2sv}}.
#' @param y vector of response data
#' @param cov1 vector of "covariate 1" values (either raw or probability 
#' transformed, whichever \code{Qcond} accepts)
#' @param cov2 vector of "covariate 2" values (either raw or probability 
#' transformed, whichever \code{Qcond} accepts)
#' @param tau_c If \code{taus=NULL}, this specifies that only quantile indices
#' larger than this will be used for regression.
#' @param K If \code{taus=NULL}, this specifies the number of equally spaced
#' quantile indices to use for regression in (tau_c, 1).
#' @param taus Keep \code{NULL} if you want to use the \code{tau_c} and
#' \code{K} option. Otherwise, specify a vector of quantile indices to
#' use for regression with.
#' @param cutoff_u Logical. If \code{TRUE}, the inputted covariates are
#' assumed to be probability transforms, and values greater than 
#' max(0.999, (1-1/(2n))) will be moved down to that maximum. As a "safety"
#' measure, this is only done if the covariates are all in the interval [0,1],
#' so if you input actual covariates that have values outside of [0,1]
#' instead of probability transforms, the truncation won't happen.
#' @return
#' Returns the objective function whose argument is the parameters (\code{cpar})
#' of the inputted \code{Qcond} function.
#' @export
Make_cnqr_obj <- function(Qcond, y, cov1, cov2,
                          tau_c=0.9, K=10, taus=NULL, cutoff_u=TRUE) {
  ## Get vector of taus, if not done already:
  if (is.null(taus)) {
    taus <- 1:K/(K+1) * (1-tau_c) + tau_c
  }
  ## If probability transforms are inputted:
  ## Ensure no u's are 1 (where some copulas are divergent). It's not realistic
  ##  anyways. Cut off at 0.999, or higher if the sample size is large
  if (cutoff_u & (all(cov1 <= 1) & all(cov2 <= 1) & all(cov1 >= 0) & all(cov2 >=0))){
    n <- length(y)
    cov1 <- sapply(cov1, function(u_) min(u_, max(0.999, (n-0.5)/n)))
    cov2 <- sapply(cov2, function(u_) min(u_, max(0.999, (n-0.5)/n)))
  }
  
  ## Output the objective function
  function(cpars){
    sum(sapply(taus, function(tau_){
      sum(rho(tau_, y - Qcond(tau_, cov1, cov2, cpars)))
    }))
  }
}