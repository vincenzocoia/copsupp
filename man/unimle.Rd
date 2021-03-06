% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/marginal.R
\name{unimle}
\alias{unimle}
\title{User-friendly univariate Maximum Likelihood Estimation}
\usage{
unimle(x, ddist = dnorm, init.val = NULL, in.range = NULL)
}
\arguments{
\item{x}{Vector of univariate sample drawn from the distribution.}

\item{ddist}{A family of densities or mass functions,
for example \code{\link{dnorm}} or \code{\link{dnbinom}}. Each parameter
should have its own named argument.}

\item{init.val}{The initial guess at parameter values for the numerical
optimizer. Could be (1) a vector, in which the order provided will correspond
to the order given in the R documentation for that distribution, or (2) a list
with names being the parameter names, or (3) \code{NULL} in which
case default parameter values are used (if there even are defaults).}

\item{in.range}{The same as the argument in \code{\link{rnlm}} -- a function
that takes the parameter values and returns \code{TRUE} if that value of
the parameter is "proper"; \code{FALSE} otherwise.}
}
\value{
A list of the parameter estimates, named according to their
corresponding parameter names
}
\description{
Finds the MLE from a univariate sample using a parametric family of
distributions built into R (or others similar to them).
}
\note{
If you want to make your own distribution (density function), just
be sure to follow R's order of arguments (\code{x} first,
then named parameters).
}
\examples{
## MLE of a Gaussian sample
unimle(rnorm(1000))
x <- rnorm(1000, 5, 100)
unimle(x)
unimle(x, init.val = list(sd = 80, mean = 10))  # Notice the order output order

## MLE of negative binomial
x <- rnbinom(1000, size = 2, mu = 10)
unimle(x, dnbinom)  # No defaults, so must specify inits
unimle(x, dnbinom, c(1, 0.3))
unimle(x, dnbinom, list(size = 1, mu = 5))
}
