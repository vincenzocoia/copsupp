% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rnlm.R
\name{rnlm}
\alias{rnlm}
\title{Restricted \code{\link{nlm}}}
\usage{
rnlm(f, p, in.range = NULL, lrg.number = 99999999, ign.error = FALSE, ...)
}
\arguments{
\item{f}{Function to minimize}

\item{p}{Starting point for the minimization}

\item{in.range}{A function that takes the
argument of \code{f} and returns TRUE if that value of the 
parameter is in the support of \code{f}, and returns FALSE if not.}

\item{lrg.number}{The large number that's put overtop of the objective
function outside of its support.}

\item{ign.error}{Sometimes \code{\link{nlm}} will throw an error (for example,
if it encounters numbers that are too big). Set this argument to \code{TRUE} to
have \code{NA} returned in the case of an error (and allow your code to keep
running afterwards).}

\item{...}{Other arguments to pass to \code{\link{nlm}}.}
}
\description{
A wrapper for \code{\link{nlm}} that puts a very large number wherever the
objective function does not exist.
}
