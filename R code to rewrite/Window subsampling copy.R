#' Subsample over a window/ball
#' 
#' Obtain the univariate subsample of responses whose covariates
#' lie in the ball around \code{x} with radius \code{r}.
#' 
#' @param xdat List of data frames as in \code{\link{ReadExisting}},
#' to be used as the covariates.
#' @param ydat List of one data frame as in \code{\link{ReadExisting}},
#' to be used as the response.
#' @param x Vector of length \code{length(xdat)}, which is the center of
#' the ball.
#' @param r Positive numeric, which is the radius of the ball.
#' @param dist Function that returns the distance of a point \code{t}
#' from 0. Default is Euclidean distance.
#' @export
#' @examples
#' spotlight(ReadExisting('totalRainfall')[2], ReadExisting('flow')[1], 2, 0.5)
spotlight <- function(xdat, ydat, x, r, dist=function(t) sqrt(sum(t^2))){
  ## Pair up.
  joint <- MatchDates.slim(xdat, ydat)
  ## Identify name of y data. Put periods where there were dashes so that it
  ##  matches the header of y in 'joint'.
  yname <- gsub("-", ".", names(ydat))
  #### Which header corresponds to y?
  yind <- which(names(joint) == yname)
  ## Subset:
  suby <- apply(joint, 1, function(row){
    ## If the x value is in the ball...
    if (dist(row[-yind]-x) <= r){
      ## ...keep the corresponding y value.
      row[yind]
    } else {
      ## ...otherwise return NA.
      NA
    }
  })
  res <- na.omit(suby)
  names(res) <- NULL
  as.vector(res) # Get rid of pesky extraneous output from na.omit.
}


#' Subsample response over a window/ball -- with non-specialized data
#' 
#' Obtain the univariate subsample of responses whose covariates
#' lie in the ball around \code{x} with radius \code{r}.
#' 
#' @param xdat Matrix of covariates (rows are observations, columns are variables)
#' @param ydat Vector of response values corresponding to the rows of \code{xdat}
#' @param x Vector of length \code{ncol(xdat)}, which is the center of
#' the ball.
#' @param r Positive numeric, which is the radius of the ball.
#' @param dist Function that returns the distance of a point \code{t}
#' from 0. Default is Euclidean distance.
#' @export
#' @examples
#' spotlight2(matrix(rnorm(3000), ncol = 3), rnorm(1000), c(0, 1, -0.5), 0.5)
spotlight2 <- function(xdat, ydat, x, r, dist=function(t) sqrt(sum(t^2))){
  ## Pair up.
  joint <- cbind(ydat, xdat)
  ## Subset:
  suby <- apply(joint, 1, function(row){
    ## If the x value is in the ball...
    if (dist(row[-1]-x) <= r){
      ## ...keep the corresponding y value.
      row[1]
    } else {
      ## ...otherwise return NA.
      NA
    }
  })
  res <- na.omit(suby)
  names(res) <- NULL
  as.vector(res) # Get rid of extraneous output from na.omit.
}