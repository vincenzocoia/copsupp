% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/augment.rvine.R
\name{augment.rvine}
\alias{augment.rvine}
\title{Add to a Column of an \code{rvine}}
\usage{
\method{augment}{rvine}(obj, a, cop, cpar, col = NULL)
}
\arguments{
\item{obj}{The object of type \code{rvine} with which to augment.}

\item{a}{Vector of vine array indices to append to the vine array column.}

\item{cop}{Vector of copula families to append to the copula matrix column.}

\item{cpar}{List of copula parameters to append to the parameter matrix column.}

\item{col}{The column of the vine to add to (i.e., the common column
of the vine array, and copula and parameter matrices). \code{NULL} (default)
if you want to add a new column to the vine.}
}
\description{
Adds extra layers to a specified column of an \code{rvine} object,
including adding a new column altogether.
}
\note{
If you're adding a new column with only one variable, there's no
need to specify the \code{cop} and \code{cpar} arguments.

Some checks are built-in to ensure that the resulting vine is indeed an
\code{rvine}, but it's not comprehensive. Would be nice to have some
function like \code{is.rvine} to check for sure.
}
\examples{
## Add to an empty vine:
rv <- rvine(matrix(nrow=0, ncol=0))
summary(rv)
rv2 <- augment(rv, a=4, col=1)
summary(rv2)

## Add to an independence vine:
rv <- rvine(matrix(4:1, ncol = 4))
summary(rv)
rv2 <- augment(rv, a=c(4,3),
               cop=c("bvtcop", "bvncop"),
               cpar=list(c(0.5, 3), -0.7),
               col=3)
summary(rv2)
summary(augment(rv2, integer(0)))

## You can't do some illegal things.
\dontrun{
augment(rv, a=1)
augment(rv, a=3:4, cop=c("frk", "frk"),
        cpar=list(4,3), col=2)
augment(rv, a=5:1, cop="frk", cpar=list(4))
}
## But you can do some:
rv2 <- augment(rv, a=c(5,5), cop="frk", cpar=list(4))
summary(rv2)
}
