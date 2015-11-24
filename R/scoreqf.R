#' Scoring functions
#'
#' Score quantile forecasts against realizations. The quantile function
#' forecasts are scored by evaluating the quantile function at multiple
#' quantile indices, instead of integrating over the functions.
#'
#' @param y Vector of realizations
#' @param yhat Forecast matrix. Each row should be a forecast with
#'  corresponding realization in \code{y}; columns are evaluations of the
#'  quantile function forecast at \code{tau}.
#' @param tau Vector of quantile indices for which quantile function is
#' evaluated.
#' @param g Increasing univariate function to transform \code{y} and \code{yhat}
#' with. Vectorized. Default is the identity.
#' @details Uses a family of proper scoring rules to score a forecast.
#' @export
scoreq <- function(y, yhat, tau, g = identity) {
  if (length(tau) != ncol(yhat))
    stop("quantile indices 'tau' and columns of yhat are not equal in number.")
  score <- 0
  for (k in 1:length(tau)) {
    score <- score + sum(rho(tau[k], g(y) - g(yhat[, k])))
  }
  score / length(y) / length(tau)
}

#' Calibration of a Forecast
#'
#' Computes the proportion of quantile forecasts that are exceeded, for each
#' quantile specified.
#'
#' @param y Vector of outcomes/response.
#' @param yhat Matrix of forecasts. Each row is a forecast corresponding to
#' \code{y}; each column corresponds to a quantile index.
#' @param tau Vector of the quantile indices forecast.
#' @param disp_plot Logical; should a P-P plot be output?
#' @return A data frame with columns "ind" for quantile indices, and "exc"
#' for corresponding proportion of outcomes that exceeded the quantile
#' of that index.
#' @export
calibration <- function(y, yhat, tau, disp_plot = TRUE) {
  if (length(tau) != ncol(yhat))
    stop("quantile indices 'tau' and columns of yhat are not equal in number.")
  res <- apply(yhat, 2, function(col){
    mean(y > col)
  })
  res <- data.frame(ind = tau, exc = res)
  if (disp_plot) {
      space <- diff(sort(tau))[1]
      p <- ggplot2::ggplot(res, ggplot2::aes(x = ind, y = exc)) +
          ggplot2::geom_abline(intercept = 1, slope = -1, linetype = "dotted") +
          ggplot2::geom_point() +
          ggplot2::scale_x_continuous("Quantile Index",
                                      limits = c(sort(tau)[1]-space, 1)) +
          ggplot2::scale_y_continuous("Proportion Exceeded",
                                      limits = c(0, max(1-(min(tau)-space), max(res$exc))))
      print(p)
  }
  res
}
