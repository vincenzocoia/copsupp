#' User-friendly univariate Maximum Likelihood Estimation
#'
#' Finds the MLE from a univariate sample using a parametric family of
#' distributions built into R (or others similar to them).
#'
#' @param x Vector of univariate sample drawn from the distribution.
#' @param ddist A family of densities or mass functions,
#' for example \code{\link{dnorm}} or \code{\link{dnbinom}}. Each parameter
#' should have its own named argument.
#' @param init.val The initial guess at parameter values for the numerical
#' optimizer. Could be (1) a vector, in which the order provided will correspond
#' to the order given in the R documentation for that distribution, or (2) a list
#' with names being the parameter names, or (3) \code{NULL} in which
#' case default parameter values are used (if there even are defaults).
#' @param in.range The same as the argument in \code{\link{rnlm}} -- a function
#' that takes the parameter values and returns \code{TRUE} if that value of
#' the parameter is "proper"; \code{FALSE} otherwise. 
#' @note If you want to make your own distribution (density function), just
#' be sure to follow R's order of arguments (\code{x} first,
#' then named parameters).
#' @return A list of the parameter estimates, named according to their
#' corresponding parameter names
#' @examples
#' ## MLE of a Gaussian sample
#' unimle(rnorm(1000))
#' x <- rnorm(1000, 5, 100)
#' unimle(x)
#' unimle(x, init.val = list(sd = 80, mean = 10))  # Notice the order output order
#'
#' ## MLE of negative binomial
#' x <- rnbinom(1000, size = 2, mu = 10)
#' unimle(x, dnbinom)  # No defaults, so must specify inits
#' unimle(x, dnbinom, c(1, 0.3))
#' unimle(x, dnbinom, list(size = 1, mu = 5))
#' @export
unimle <- function(x, ddist = dnorm, init.val = NULL, in.range = NULL) {
  ## Get names of parameters
  if (is.null(init.val)) {
    ## init.val not specified. Take the defaults.
    ## Note: Assume (besides x) the only non-parameter option is possibly "log"
    #### Get arguments
    arg <- formals(ddist)
    arg <- as.list(arg) # Change from "Dotted" list to list.
    #### Get rid of log and "x"
    arg <- arg[-1]
    arg$log <- NULL
    #### Find out what arguments remain that are specified.
    keep <- sapply(arg, is.vector)
    #### Use those defaults as the initial values
    init.val <- arg[keep]
    #### Capture parameter names
    parnames <- names(init.val)
  }
  if (is.vector(init.val)) {
    ## Initial values are inputted as a vector. Take the first non-x
    ##  arguments to be the parameters.
    parnames <- formalArgs(ddist)[1+1:length(init.val)]
  }
  if (is.list(init.val)) {
    ## User entered list of parameters, supposedly named.
    parnames <- names(init.val)
    if (length(parnames) != length(init.val)){
      stop("If you enter a list for init.val, it must be named by the parameter names.")
    }
  }

  ## Set up a mapping from a vector of parameters to arguments of ddist.
  vec2arg <- function(vec){
    res <- as.list(vec)
    names(res) <- parnames
    c(list(x = x), res)
  }
  ## Finally, set up neg. log likelihood.
  nllh <- function(theta){
    arg <- vec2arg(theta)
    -sum(log(do.call(ddist, arg)))
  }
  ## Minimize:
  fit <- rnlm(nllh, c(init.val, recursive = T), in.range = in.range)
  vec2arg(fit$estimate)[-1]
}


#' Estimate Marginal distribution
#'
#' Estimates the marginal distribution of a univariate sample, either
#' by MLE (if a parametric family is specified) or empirical distribution (edf,
#' with kernel density pdf)
#' otherwise.
#'
#' @param x Vector representing the univariate sample drawn from the
#' distribution to be estimated.
#' @param dist.name String; name of the parametric family of distributions, if
#' you don't want the edf. Use R's naming conventions
#' for distributions ahead of the d/p/q/r. For example, the Normal
#' distribution would be \code{"norm"}.
#' @param ecdf.split Numeric; real number for which to fit the specified model
#' for data above such number, and empirical distribution to be fit below.
#' @param init.val If a parametric distribution is specified (through
#' \code{dist.name}), this should be the initial value specification for the
#' parameters. Could be (1) a vector, in which the order provided will correspond
#' to the order given in the R documentation for that distribution, or (2) a list
#' with names being the parameter names, or (3) remain \code{NULL} in which
#' case default parameter values are used (if there even are defaults).
#' @param soften.cdf Logical; should the cdf be modified so that it never returns
#' the boundary cases 0 or 1? (uses 0.5/n or 1-0.5/n instead,
#' where n is amount of data)
#' @param eqf.type If the empirical quantile function is being estimated, this
#' argument specifies the type of quantile algorithm to use, as in the
#' \code{type} argument of the \code{\link{quantile}} function. In particular,
#' an integer between 1 and 9.
#' @param in.range The same as the argument in \code{\link{rnlm}} -- a function
#' that takes the parameter values and returns \code{TRUE} if that value of
#' the parameter is "proper"; \code{FALSE} otherwise. 
#' @details Some noteworthy options for \code{eqf.type} are: 1 = inverse
#' empirical cdf (default of \code{marginal}); 7 = default of
#' \code{\link{quantile}}, a continuous function.
#' @return If \code{return.param = FALSE} (default), output is a list of
#' length three:
#'
#' \code{$cdf}: The estimated distribution function
#'
#' \code{$qf}: The estimated quantile function (inverse of the cdf estimate)
#'
#' \code{$pdf}: The estimated density function or mass function. Gives a kernel
#' density estimate if using the empirical distribution.
#' @export
marginal <- function(x,
                     dist.name = NULL,
                     ecdf.split = NULL,
                     init.val = NULL,
                     in.range = NULL,
                     soften.cdf = TRUE,
                     eqf.type = 1){
  numorig <- length(x)
  x <- na.omit(x)
  n <- length(x)
  if (n != numorig)
      warning(paste("Function 'marginal()' is ignoring the", numorig - n, "NA values."))
  ## Case 1: no model given. Use empirical distribution.
  if (is.null(dist.name)){
    ## Get edf if no distribution family is given.
    ## edf
    cdf <- ecdf(x)
    attributes(cdf) <- NULL
    ## quantile function (it's vectorized over tau)
    qf <- function(tau){
      res <- quantile(x, tau, type = eqf.type)
      attributes(res) <- NULL
      res
    }
#     ## mass function... for completeness only.
#     tab <- table(x)
#     nams <- names(tab)
#     nx <- length(x)
#     pdf <- function(q){
#       mass <- sapply(q, function(q_){
#         res <- which(nams == as.character(q_))
#         if (length(res) == 0){
#           res <- 0
#         } else {
#           res <- tab[res] / nx
#         }
#         res
#       })
#       names(mass) <- NULL
#       mass
#     }
    ## Use kernel density instead of mass function, since we're probably working
    ##   with continuous variables.
    #### Density estimate. Use default (only return points, not a function)
    dens <- density(x)
    pdf <- approxfun(dens$x, dens$y)
  } else {  ## Case 2: Model is given
    ## Get distribution names
    ddist <- get(paste0("d", dist.name))
    qdist <- get(paste0("q", dist.name))
    pdist <- get(paste0("p", dist.name))
    ## Just to be safe, get the name of the first argument of those functions:
    d1 <- formalArgs(ddist)[1]
    q1 <- formalArgs(qdist)[1]
    p1 <- formalArgs(pdist)[1]
    ## Case 2a: Use the model for all the data:
    if (is.null(ecdf.split)) {
        ## Get estimate
        thetahat <- unimle(x, ddist, init.val = init.val, in.range = in.range)
        ## Get distribution-related functions:
        cdf <- function(q){
            arg <- c(list(q), thetahat)
            names(arg)[1] <- p1
            do.call(pdist, arg)
        }
        pdf <- function(x){
            arg <- c(list(x), thetahat)
            names(arg)[1] <- d1
            do.call(ddist, arg)
        }
        qf <- function(tau){
            arg <- c(list(tau), thetahat)
            names(arg)[1] <- q1
            do.call(qdist, arg)
        }
    } else {  ## Case 2b: This model only applies above some threshold.
        ## Split the data
        xlow <- subset(x, x < ecdf.split)
        xhigh <- subset(x, x >= ecdf.split)
        ## Estimate probability of being low or high:
        plow <- length(xlow) / length(x)
        phigh <- 1 - plow
        ## Get lower and upper fits:
        lowfit <- marginal(xlow, 
                           eqf.type = eqf.type, 
                           soften.cdf = soften.cdf, 
                           in.range = in.range)
        highfit <- marginal(xhigh, 
                            dist.name = dist.name, 
                            init.val = init.val, 
                            soften.cdf = soften.cdf,
                            in.range = in.range)
        ## Combine.
        cdf <- function(x) {
            sapply(x, function(x_) {
                if (x_ < ecdf.split) {
                    lowfit$cdf(x_) * plow
                } else {
                    plow + highfit$cdf(x_) * phigh
                }
            })
        }
        pdf <- function(x) {
            lowfit$pdf(x) * plow +
                highfit$pdf(x) * phigh
        }
        qf <- function(u) {
            sapply(u, function(u_) {
                if (u_ < plow){
                    lowfit$qf(u_ / plow)
                } else {
                    highfit$qf((u_ - plow) / phigh)
                }
            })
        }
    }

  }
  ## Lastly, if softened edges are wanted, make them:
  if (soften.cdf) {
      hardcdf <- cdf
      ## Remake cdf:
      cdf <- function(x) {
          ans <- hardcdf(x)
          sapply(ans, function(ans_) max(min(1-0.5/n, ans_), 0.5/n))
      }
  }
  list(cdf = cdf, qf = qf, pdf = pdf)
}
