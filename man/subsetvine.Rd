% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/subset.rvine.R
\name{subset.rvine}
\alias{subset.rvine}
\alias{subsetvarray}
\title{Subset a Regular Vine}
\usage{
\method{subset}{rvine}(rv, select, justcheck = FALSE)

subsetvarray(G, select, justcheck = FALSE)
}
\arguments{
\item{rv}{A regular vine object.}

\item{select}{Vector of variables to subset,
if possible. The order of the variables does not matter.}

\item{justcheck}{Logical; should this function only check whether or not
the subset exists? \code{TRUE} if so.}

\item{G}{A vine array matrix.}
}
\value{
If \code{justcheck} is \code{TRUE}, returns \code{TRUE} if the requested
subset exists, and \code{FALSE} if not.

If \code{justcheck} is \code{FALSE}, returns
a vine of the subsetted variables, with variables
ordered according to their order in \code{G}; or \code{NULL} if
the subset does not form a vine.
}
\description{
Find the vine connecting a subset of variables from a bigger vine, if
it exists. \code{subsetvarray} only finds the subset of a vine array.
}
\details{
Just a technicality:
by saying a subset "doesn't have an existing vine", I mean that
a vine can't be formed using nodes and edges from
the original -- not that the
joint distribution of the selected variables can't be created from a vine
(so as to say, for example, that the simplifying assumption of vines
doesn't hold for this distribution).
}
\examples{
## Setup a vine.
G <- AtoG(CopulaModel::Dvinearray(5))
subsetvarray(G, c(2, 4, 3))

copmat <- makeuppertri(c("gum", "mtcj", "gal", "joe",
                         "frk", "gum", "bb7",
                         "bb1", "indepcop",
                         "bb8"), 4, 5, "")
cparmat <- makeuppertri.list(c(3, 2.5, 2, 1.5,
                               1, 1.3, 2, 2,
                               3, 4,
                               5, 0.5),
                               len = c(1,1,1,1,1,1,2,2,0,2),
                               4, 5)
(rv <- rvine(G, copmat, cparmat))

## Subset some variables.
subset(rv, c(2, 4, 3))
subset(rv, 5)
subset(rv, integer(0))

## This subset won't work:
subset(rv, c(4, 1), justcheck = TRUE)
## But it will in a 0-truncated vine:
subset(trunc(rv, 0), c(4, 1), justcheck = TRUE)
subset(trunc(rv, 0), c(4, 1))

## Select variables not present?
subset(rv, c(2, 4, 17))
}
