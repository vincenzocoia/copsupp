#' Quantile Model Predictions
#'
#' @param appendy A fitted object with \code{\link{fitlayer_cnqr}}.
#' @param dat Either a data matrix to make predictions on, or one of the
#' two strings: \code{"val"} will make predictions on the validation set, and
#' \code{"tr"} will make predictions on the training set.
#' @param tau Vector of quantile indices to make predictions at.
#' @export
predict_appendy <- function(appendy, dat = "val", tau = space_taus(10)) {
    ## Get sequential conditional preidtors.
    ucond <- NULL
    if (is.character(dat)) {
        if (dat == "val") {
            ucond <- appendy$ucondval
        }
        if (dat == "tr") {
            ucond <- appendy$ucondtr
        }
    } else {
        v <- appendy$xord
        rv <- appendy$basevine
        ucond <- pcondseq(dat, v, rv)
    }
    if (is.null(ucond))
        stop("'dat' input is inadequate.")
    ## Forecast.
    cops <- appendy$cops
    cpars <- appendy$cpars
    QY <- appendy$QY
    QYgX(tauset, ucond, cops = cops, cpars = cpars, QYfitted = QY)
}
