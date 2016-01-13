#' Create an object of class `rvine`
#'
#' @param G G-Vine array
#' @param copmat Upper triangular \code{(nrow(G)-1) x ncol(G)} matrix of copula model names.
#' Or, a vector of length 1 if the same copula model applies to all
#' edges.
#' @param cparmat Upper triangular \code{(nrow(G)-1) x ncol(G)} matrix of
#' parameters taken by the copulas in \code{copmat}. Or, a single vector
#' of the parameter to be applied to the entire vine.
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
#' G <- AtoG(CopulaModel::Dvinearray(4))
#' (rv <- rvine(G, "frk", 2))
#' lapply(rv, identity)
#' (rv <- rvine(G, "bvtcop", c(0.4, 5)))
#'
#' G <- AtoG(CopulaModel::Dvinearray(4))
#' copmat <- makevinemat("gum",
#'                       c("bvtcop", "frk"),
#'                       c("mtcj", "frk", "indepcop"),
#'                       c("bvncop", "joe", "mtcj", "frk"), zerocol = TRUE)
#' cparmat <- makevinemat(3.1,
#'                        list(c(0.5, 4), 2.3),
#'                        list(4.2, 3.5, numeric(0)),
#'                        c(0.5, 2.2, 2.5, 1.6), zerocol = TRUE)
#' rv <- rvine(G, copmat, cparmat)
#' @export
rvine <- function(G, copmat, cparmat) {
    ## First -- deal with the "trivial case" of G.
    d <- ncol(G)
    if (d == 0) {
        return(structure(list(G=G, copmat=NA, cparmat=NA), class = "rvine"))
    }
    if (d == 1) {
        return(structure(list(G=G,
                              copmat=matrix(nrow=0, ncol=1),
                              cparmat=matrix(nrow=0, ncol=1)),
                         class = "rvine"))
    }
    ntrunc <- nrow(G) - 1
    if (ntrunc == 0) {
        return(structure(list(G=G,
                              copmat=matrix(nrow=0, ncol=d),
                              cparmat=matrix(nrow=0, ncol=d)),
                         class = "rvine"))
    }
    ## Next, construct the copula matrix if not already done.
    if (!is.matrix(copmat)) {
        copmat <- rep(copmat, ntrunc*d - choose(ntrunc+1, 2))
        copmat <- makeuppertri(copmat, ntrunc, d, blanks = "")
    }
    ## Copula parameter matrix:
    if (is.matrix(cparmat)) if(!is.list(cparmat[1,1])) {
        cparvec <- c(t(cparmat)[lower.tri(t(cparmat))], recursive = TRUE)
        cparmat <- makeuppertri.list(cparvec, len=rep(1,length(cparvec)),
                                     nrow = nrow(cparmat), ncol = ncol(cparmat))
    }
    if (!is.matrix(cparmat)) {
        cparmat <- rep(cparmat, ntrunc*d - choose(ntrunc+1, 2))
        len <- rep(length(cparmat), ntrunc*d - choose(ntrunc+1, 2))
        cparmat <- makeuppertri.list(cparmat, len, nrow=ntrunc, ncol=d)
    }
    ## Output
    res <- list(G = G,
                copmat = copmat,
                cparmat = cparmat)
    class(res) <- "rvine"
    res
}

#' @export
print.rvine <- function(rv) {
    d <- ncol(rv$G)
    if (d == 0) return(cat("Empty vine: no variables."))
    v <- rv$G[1, ]
    ntrunc <- nrow(rv$G) - 1
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
    invisible()
}

#' @export
is.rvine <- function(rv) {
    inherits(rv, "rvine")
}
