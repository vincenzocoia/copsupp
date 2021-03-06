% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/makeuppertri.R
\name{makeuppertri}
\alias{makeuppertri}
\alias{makeuppertri.list}
\title{Make an upper triangular matrix}
\usage{
makeuppertri(
  entries,
  row,
  col,
  blanks = 0,
  byRow = TRUE,
  incDiag = FALSE,
  nrow = NULL,
  ncol = NULL
)

makeuppertri.list(
  entries,
  len,
  row,
  col,
  blanks = list(),
  byRow = TRUE,
  incDiag = FALSE,
  nrow = NULL,
  ncol = NULL
)
}
\arguments{
\item{entries}{Vector of entries for the upper triangular part of the matrix.}

\item{row}{Number of rows in the matrix}

\item{col}{Number of columns in the matrix}

\item{blanks}{What should go in the non-diagonal entries?}

\item{byRow}{Logical; read the \code{entries} by row (as if you're reading)
(TRUE), or read in vertically (FALSE).}

\item{incDiag}{Should the entries go on the diagonal too? \code{TRUE} if
so, \code{FALSE} if not.}

\item{nrow, ncol}{When I first made this function, I used \code{row} and
\code{col}, but I should have called them \code{nrow} and \code{ncol}.
Hence this addition.}

\item{len}{Vector of positive integers which specify the lengths of the
individual vectors that are pooled in \code{entries}.}
}
\description{
Upper triangular matrix making made easy.
}
\note{
Use \code{makeuppertri} to make a matrix. If you want entries to be
vectors (which would have to be an array with list entries), use
\code{makeuppertri.list}.
}
\examples{
## Square matrices
makeuppertri(1:choose(5,2), 5, 5)
makeuppertri(1:choose(6,2), 5, 5, incDiag = TRUE)

## Not square.
makeuppertri(1:9, row = 2, col = 5, incDiag = TRUE)
makeuppertri(1:3, row = 5, col = 3)
## List entries
(M <- makeuppertri.list(1:12, c(1, 10, 1), 3, 3))
M[1, 2]   # A list of length one.
M[1, 3]   # Another list of length one.
(M <- makeuppertri.list(1:3, rep(1,3), 3, 3))
M[1, 2]   # Still a list of length one.
}
