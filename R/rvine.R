#' Create an object of class `rvine`
#'
#' @param A Vine array (matrix), possibly truncated.
#' @param copmat Upper triangular \code{(nrow(A)-1) x ncol(A)} matrix of copula model names.
#' Or, a vector of copula model names for the edges in \code{A}, listed in
#' "reading order". Use a vector of length 1 if the same copula model applies to all
#' edges.
#' @param cparmat Upper triangular \code{(nrow(A)-1) x ncol(A)} matrix of
#' parameters taken by the copulas in \code{copmat}.
#' @param marg List of \code{ncol(A)} vectorized marginal distribution functions
#' with entries corresponding to the columns in \code{A}. Or, a single vectorized
#' (distribution) function that applies to each variable.
#' @return An object of class `rvine`, which is a named list of the arguments.
#' @examples
#' ## Empty vine:
#' rvine(matrix(integer(0)))
#'
#' ## Independence vine:
#' (rv <- rvine(matrix(4:1, ncol = 4)))
#' ## Take a look at each component of the vine:
#' lapply(rv, identity)
#'
#' A <- CopulaModel::Dvinearray(4)
#' rvine(A, "frk", makeuppertri(2, nrow=3, ncol=4), pexp)
#' rvine(A, rep(c("frk", "gum"), 3), makeuppertri(2, nrow=3, ncol=4))
#'
#' copmat <- makeuppertri(c("gum", "bvtcop", "mtcj",
#'                          "frk", "indepcop",
#'                          "frk"), 3, 4, blanks = "")
#' rvine(A, copmat)
#'
#' cparmat <- makeuppertri.list(c(2, 0.5, 4, 2,
#'                                1,
#'                                1), len = c(1,2,1,1,0,1), nrow = 3, ncol = 4)
#' rvine(A, copmat, cparmat)
#' @export
rvine <- function(A, copmat = NULL, cparmat = NULL, marg = identity) {
    ## First -- deal with the "trivial case" of A.
    d <- ncol(A)
    if (d == 0) {
        return(structure(list(A=A, copmat=NA, cparmat=NA, marg=NA), class = "rvine"))
    }
    if (d == 1) {
        return(structure(list(A=A, copmat=matrix(nrow=0, ncol=1),
                         cparmat=matrix(nrow=0, ncol=1), marg = marg),
                         class = "rvine"))
    }
    ## Next, construct the copula matrix if not already done.
    ntrunc <- nrow(A) - 1
    if (ntrunc == 0){
        copmat <- matrix(nrow = 0, ncol = d)
        cparmat <- matrix(nrow = 0, ncol = d)
    } else {
        if (!is.matrix(copmat) & !is.null(copmat)) {
            if (length(copmat) == 1) {
                copmat <- rep(copmat, ntrunc*d - choose(ntrunc+1, 2))
            }
            copmat <- makeuppertri(copmat, ntrunc, d, blanks = "")
        }
    }
    ## Copula parameter matrix:
    if (is.null(copmat) & !is.null(cparmat)) {
        warning("Parameter matrix cannot be specified until copula matrix is.")
        cparmat <- NULL
    }
    if (is.matrix(cparmat)) if (!is.list(cparmat[1,1])) {
        cparvec <- c(t(cparmat)[lower.tri(t(cparmat))], recursive = TRUE)
        cparmat <- makeuppertri.list(cparvec, len=rep(1,length(cparvec)),
                                     nrow = nrow(cparmat), ncol = ncol(cparmat))
    }
    ## Make vector of cdfs if it's a single value:
    if (length(marg) == 1) marg <- rep(list(marg), d)
    ## Output
    res <- list(A = A,
                copmat = copmat,
                cparmat = cparmat,
                marg = marg)
    class(res) <- "rvine"
    res
}

#' @export
print.rvine <- function(rv) {
    d <- ncol(rv$A)
    if (d == 0) return(cat("Empty vine: no variables."))
    v <- vars(rv)
    ntrunc <- nrow(rv$A) - 1
    if (ntrunc == 0) {
        trunctext <- "Independent"
    } else {
        trunctext <- paste0(ntrunc, "-truncated")
    }
    cat(paste0(trunctext, " vine with variables ", paste(v, collapse = ", "), ".\n"))
    if (ntrunc > 0) {
        if (!is.null(rv$copmat) & is.null(rv$cparmat))
            cat("\nUnspecified parameters.")
        if (is.null(rv$copmat))
            cat("\nUnspecified copula models and parameters.")
    }
    if (identical(rv$marg, rep(list(identity), d)))
        cat("\nUniform margins.")
    invisible()
}
