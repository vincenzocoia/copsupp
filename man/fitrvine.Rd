% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fitrvine.R
\name{fitrvine}
\alias{fitrvine}
\title{Fit a Vine Model}
\usage{
fitrvine(dat, layer = layeropts(1:ncol(dat)), basevine = NULL, ...)
}
\arguments{
\item{dat}{Matrix of data having Uniform marginal distributions;
columns represent variables, and rows observations.}

\item{...}{Other arguments to pass to \code{VineCopula::RVineCopSelect}.}

\item{rv}{Object of class "rvine", representing a pre-specification of the
vine to fit. Or, \code{NULL} for fully unspecified. See details.}

\item{var}{Vector of integers specifying the column numbers of \code{dat}
to fit a model to. Default is all variables.}

\item{ntrunc}{Integer, either \code{1, 2, ...,ncol(dat)-1},
of the truncation level of the vine to be fit.}

\item{families}{A vector of copula family names to try
fitting (will also consider their rotations/reflections). Limited to
those families available in \code{VineCopula} package, listed in
\code{\link{BiCopSelect}}.}
}
\value{
A "fitrvine" object, which has class \code{c("fitrvine", "rvine")},
which is a named list of the following:

\itemize{
     \item \code{$A}: Vine array, truncated to \code{ntrunc}.
     \item \code{$copmat}: \code{ntrunc x ncol(A)} upper-triangular
     matrix of copula model names.
     \item \code{$cparmat}: \code{ntrunc x ncol(A)} upper-triangular
     matrix of copula parameters. Each entry is a list of length one containing
     the vector of copula parameters for that copula family.
     \item \code{$dat}: The inputted data matrix, \code{dat}.
     \item \code{$aic}: The AIC of the fitted model.
     \item \code{$bic}: The BIC of the fitted model.
     \item \code{$nllh}: The negative log likelihood of the fitted model.
     \item \code{$covmat}: Covariance matrix of the fitted parameters (in
     reading-order of the parameters, i.e.
     \code{c(t(cparmat)[lower.tri(t(cparmat))], recursive=TRUE)}).
}
}
\description{
Fits a joint distribution for the data using a vine copula. The vine array
is first chosen using the minimum spanning tree algorithm using
the function \code{gausstrvine.mst} in the \code{CopulaModel}
package, then the pairwise
copula models are chosen and fit using \code{RVineCopSelect} in the
\code{VineCopula} package.
}
\details{
If you want to specify parts of the vine, then specify them in an "rvine" object
(see \code{\link{rvine}}):

\enumerate{
     \item Your first option is to specify the vine array \code{A}.
     \item If the vine array is specified, then you can specify some or all
     of the copula families by putting them in \code{copmat}. Leave
     unspecified edges as \code{NA}.
     \item If there are copula families specified, you can specify parameters
     for those families by putting the parameters in \code{cparmat}.
     Unspecified parameters should be \code{NA}.
}

For parts of the vine that are unspecified, you have some fitting options:

\itemize{
     \item If you didn't specify a vine array \code{A}, you can select which
     variables (column numbers of \code{dat}) you'd like to fit through the
     \code{var} argument. You can also select the truncation level of
     the vine array through the argument \code{ntrunc}.
     \item If you left some copula families unspecified, you can indicate
     the candidate families in the \code{families} argument.
}
}
\examples{
## Get some simulated data:
set.seed(152)
ntrunc <- 2
d <- 4
A0 <- truncvarray(CopulaModel::Dvinearray(d), ntrunc)
copmat0 <- makeuppertri("frk", ntrunc, d, "")
cparmat0 <- makeuppertri(3, ntrunc, d)
dat <- fvinesim(100, A0, copmat0, cparmat0)

## Fit a model to the data:
fit.rvine(dat, ntrunc=ntrunc)
fit.rvine(dat, c(4, 2, 3))
}
