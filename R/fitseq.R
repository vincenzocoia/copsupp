#' Choose sequence of variables to link
#'
#' Chooses the sequence of variables that link up with another variable
#' in the order of highest to lowest partial correlation, as approximated by
#' lm(). Intended to internal use. Use \code{fitseq_rvine} if you want
#' the order to be chosen so that, when attached to a base rvine, you end
#' up with an rvine still (NOT FUNCTIONAL).
#'
#' @param dat Matrix of data with Unif(0,1) margins
#' @param var Integer; the variable you want to link.
#' @param linkwith Integer vector; the variables you want to link up \code{var}
#' with.
#' @param a Vector of a pre-specified vine array column, if you know
#' some of it. Should start with \code{var}, and have \code{NA} for blank spots.
#' @details
#' All variables in \code{linkwith} will be linked, despite what's in \code{a}.
#' @return Vector of the variable sequence, with \code{var} first.
#' @examples
#' set.seed(123)
#' dat <- matrix(runif(50), ncol=5)
#' fitseq(dat, 1, 3:5)
#' fitseq(dat, 1, 2:5, a=c(1,NA,4))
#' @rdname fitseq
#' @export
fitseq <- function(dat, var, linkwith = (1:ncol(dat))[-var], a = NULL) {
    if (!is.null(a)) linkwith <- setdiff(linkwith, a[-1])
    if (is.null(a)) a <- c(var, rep(NA, length(linkwith)))
    if (sum(is.na(a)) < length(linkwith)) a <- c(a, rep(NA, length(linkwith) - sum(is.na(a))))
    for (i in 2:length(a)) if (is.na(a[i])) {
        condset <- a[1+seq_len(i-2)]
        candidatevars <- setdiff(linkwith, condset)
        y <- qnorm(dat[, var])
        x <- qnorm(dat[, candidatevars])
        if (!is.matrix(x)) x <- matrix(x, ncol = 1)
        cond <- qnorm(dat[, condset])
        if (length(condset) != 0) {
            y <- lm(y ~ cond)$residuals
            names(y) <- NULL
            x <- apply(x, 2, function(x_) {
                res <- lm(x_ ~ cond)$residuals
                names(res) <- NULL
                res
            })
        }
        cors <- cor(x, y)
        a[i] <- candidatevars[which(cors == max(cors))]
    }
    a
}

#' @param A Vine array matrix containing variables in \code{linkwith} and
#' \code{a[-1]} (to ensure the order is chosen to preserve the rvine status).
#' @rdname fitseq
#' @export
fitseq_rvine <- function(dat, var, linkwith = (1:ncol(dat))[-var], a=NULL, A) {
#     if (!is.null(a)) linkwith <- setdiff(linkwith, a[-1])
#     if (is.null(a)) a <- c(var, rep(NA, length(linkwith)))
#     if (sum(is.na(a)) < length(linkwith)) a <- c(a, rep(NA, length(linkwith) - sum(is.na(a))))
#     for (i in 2:length(a)) if (is.na(a[i])) {
#         condset <- a[1+seq_len(i-2)]
#         candidatevars <- setdiff(linkwith, condset)
#         y <- qnorm(dat[, var])
#         x <- qnorm(dat[, candidatevars])
#         if (!is.matrix(x)) x <- matrix(x, ncol = 1)
#         cond <- qnorm(dat[, condset])
#         if (length(condset) != 0) {
#             y <- lm(y ~ cond)$residuals
#             names(y) <- NULL
#             x <- apply(x, 2, function(x_) {
#                 res <- lm(x_ ~ cond)$residuals
#                 names(res) <- NULL
#                 res
#             })
#         }
#         cors <- cor(x, y)
#         a[i] <- candidatevars[which(cors == max(cors))]
#     }
#     a
}
