% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fitlayer.R
\name{fitlayer}
\alias{fitlayer}
\title{Fit a new layer}
\usage{
fitlayer(
  dat,
  basevine,
  edges,
  cops = NULL,
  cpars = NULL,
  families = c("indepcop", "bvncop", "bvtcop", "mtcj", "gum", "frk", "joe", "bb1",
    "bb7", "bb8")
)
}
\arguments{
\item{dat}{Data matrix with Uniform margins.}

\item{basevine}{Object of type "rvine" of the already-fit base vine for which
the new layer is to be applied.}

\item{edges}{Vector -- new column of vine array (with node appearing first)}

\item{cops}{Vector or list of pre-specified
copula families for each edge. Put \code{NA}
to leave the edge unspecified. \code{NULL} for
fully unspecified. You're allowed to put more than one family
as candidates.}

\item{cpars}{Pre-specified copula parameters corresponding to some of the
specified copulas in \code{cops}. Put \code{NA} in place of parameters to
leave them unspecified. \code{NULL} for fully unspecified.}

\item{families}{Vector of candidate copula family names for those that are
\code{NA} or \code{NULL}.}
}
\value{
List of fitted \code{$cops} and \code{$cpars}.
}
\description{
Choose and fit copula models on a new layer. The edge ("array column") must be
pre-specified. Intended for internal use.
}
\details{
Edges are fit so that edges[1] is the "V" variable. So, copulas
are fit to (edges[2], edges[1]), then (edges[3],
edges[1]) | edges[2], etc. That's because when computing edges[1]|others,
"pcond" can be used instead of "pcond12".
}
\note{
Expecting smart input. So, ensure that \code{edges} has length at
least 2, and that edges[-1] are variables in \code{basevine}, and that
\code{cpars} are only specified when there's only one copula family to
choose from.
}
