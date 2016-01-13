#' Convert \code{"rvine"} to \code{"fitrvine"}
#'
#' Converts a fully specified regular vine into a fitted regular vine after
#' specifying data. Since the original regular vine is a fully specified
#' distribution, there are no fitted parameters.
#'
#' @param dat Matrix of data with Unif(0,1) margins.
#' @param basevine Object of type 'rvine' (completely specified).
#' @return Object of type 'fitrvine'.
#' @export
rvine2fitrvine <- function(dat, basevine) {
    nvar <- length(vars(basevine))
    basevine$dat <- dat
    basevine$estimated <- basevine$cparmat
    if (nvar == 0) {
        basevine$nllh <- NA
        basevine$aic <- NA
        basevine$bic <- NA
        class(basevine) <- c("fitrvine", "rvine")
        return(basevine)
    }
    if (nvar == 1) {
        ## AIC/BIC are 0 because the density is always 1.
        basevine$nllh <- 0
        basevine$aic <- 0
        basevine$bic <- 0
        class(basevine) <- c("fitrvine", "rvine")
        return(basevine)
    }
    ## Put FALSE in place of parameters in "estimated" matrix.
    basevine$estimated <- apply(basevine$estimated, 1:2, function(l) {
        if (is.null(l[[1]])) NULL else rep(FALSE, length(l[[1]]))
    })
#     basevine$estimated <- apply(basevine$estimated, 1:2, function(l) {
#         if (is.null(l[[1]])) {
#             NULL
#         } else {
#             sapply(l[[1]], function(v_) if (is.na(v_)) TRUE else FALSE)
#         }
#     })
    ## Get likelihood
    nllh <- -sum(logdrvine(dat, basevine))
    basevine$nllh <- nllh
    basevine$aic <- 2 * nllh
    basevine$bic <- 2 * nllh
    class(basevine) <- c("fitrvine", "rvine")
    basevine
}
