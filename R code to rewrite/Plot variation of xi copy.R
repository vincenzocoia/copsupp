#' Plotting data of conditional EVI using moving window
#' 
#' For univariate x and y, a moving-window is applied to estimate
#' the conditional EVI across the range of x values. Returns a data
#' frame of the x value and window width, along with the corresponding 
#' conditional EVI estimates and estimated standard
#' errors. 
#' 
#' @param xdat List of one data frame from \code{\link{ReadExisting}} to be used
#' as the x variable.
#' @param ydat List of one data frame from \code{\link{ReadExisting}} to be used
#' as the y variable.
#' @param w Positive numeric; specifies the fixed width of the window to use. If
#' \code{NULL}, chooses the window width that ensures the overlaying subsample
#' is of size \code{ensure_n}.
#' @param ensure_n Positive integer; the number of observations that should overlay
#' the window. Only used if \code{w=NULL} (i.e. non-fixed window
#' width)
#' @param Nx Positive integer; number of points along the x-axis to estimate
#' the EVI at.
#' @param estim Univariate estimator of the EVI, such as \code{\link{evi.Hill}}.
#' @param ... arguments to pass to the \code{estim} function.
#' @export
xiplotdat <- function(xdat, ydat, w = NULL, ensure_n = 100, 
                          Nx = 100, estim = evi.Hill, ...){
  library(ggplot2)
  library(plyr)
  ## Pair up the data (if anything, necessary for subsampling so that
  ##  the dates match).
  paired <- MatchDates.slim(xdat, ydat)
  xname <- gsub("-", ".", names(xdat))
  #### Make sure the asked-for subsample size is no more than the amount of
  ####  available data.
  ensure_n <- min(ensure_n, nrow(paired))
  #### Don't bother continuing if there's no data in the first place.
  if (nrow(paired) == 0)
    return(data.frame(x=numeric(0), window=numeric(0), estimate=numeric(0), SE=numeric(0)))
  ## Get set of xvalues for which to estimate EVI index at.
  xsample <- sort(paired[[xname]])
  xwidth <- tail(xsample, 1) - head(xsample, 1)
  xset <- xsample[1] + (1:Nx) * xwidth / (Nx+1)
  #### Put this in the to-be plotting data frame.
  plotdat <- data.frame(x = xset)
  ## Get corresponding window widths, if they're not specified (and thus fixed)
  if (is.null(w)){
    for (x in xset){
      w <- c(w, sort(abs(xsample - x))[ensure_n] * 2)
    }
  } else {
    w <- rep(w, Nx)
  }
  #### Put window widths in plotting data frame.
  plotdat <- cbind(plotdat, window = w)
  ## Compute estimates and SE's over each window.
  plotdat <- ddply(plotdat, names(plotdat), function(dfrow){
    ## Get subsample of y's using ball function
    suby <- spotlight(xdat, ydat, x=unique(dfrow$x), r=unique(dfrow$window)/2)
    ## Get estimate EVI and get SE
    if (!exists("k")){
      k <- floor(length(suby) * 0.3)
      message(paste("k not given for Hill estimator. Using", k))
    }
    ## Need to change k if using Pickands estimator, since it uses 4k+1, not k+1.
    est <- estim(suby, k)
    if (all(is.na(est)))
      est <- data.frame(estimate=NA, SE=NA)
    est
  })
  ## Output plotting data
  plotdat
}