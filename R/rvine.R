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
#' rvine(matrix(nrow=0, ncol=0))
#'
#' ## Independence vine:
#' (rv <- rvine(matrix(4:1, ncol = 4)))
#' ## Take a look at each component of the vine:
#' summary(rv)
#'
#' G <- AtoG(CopulaModel::Dvinearray(4))
#' (rv <- rvine(G, "frk", 2))
#' summary(rv)
#' (rv <- rvine(G, "bvtcop", c(0.4, 5)))
#'
#' G <- AtoG(CopulaModel::Dvinearray(5))
#' copmat <- makevinemat("gum",
#'                       c("bvtcop", "frk"),
#'                       c("mtcj", "frk", "indepcop"),
#'                       c("bvncop", "joe", "mtcj", "frk"),
#'                       zerocol = TRUE)
#' cparmat <- makevinemat(3.1,
#'                        list(c(0.5, 4), 2.3),
#'                        list(4.2, 3.5, numeric(0)),
#'                        c(0.5, 2.2, 2.5, 1.6), zerocol = TRUE)
#' rv <- rvine(G, copmat, cparmat)
#' summary(rv)
#' @export
rvine <- function(G, copmat, cparmat) {
    ## First -- deal with the "trivial case" of G.
    d <- ncol(G)
    if (d == 0) {
        return(structure(list(G=matrix(ncol=0, nrow=1), copmat=NA, cparmat=NA),
                         class = "rvine"))
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
        len <- rep(length(cparmat), ntrunc*d - choose(ntrunc+1, 2))
        cparmat <- rep(cparmat, ntrunc*d - choose(ntrunc+1, 2))
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

#' @export
summary.rvine <- function(rv) {
    summat <- combine_copmat(rv)
    cat("Vine Array:\n")
    print(rv$G)
    cat("\nCopulas:\n")
    print(summat)
}

#' Combine Copula and Parameter Matrices
#'
#' Combines a copula matrix and a parameter matrix, so that specific
#' copulas are indicated.
#'
#' @param rv Object of type "rvine"
#' @param copmat Matrix of copula families, as in an "rvine" object.
#' @param cparmat Matrix of copula parameters, as in an "rvine" object.
#' @param digits Number of significant digits to round parameter values to.
#' @note \code{copmat} and \code{cparmat} are not needed if \code{rv} is specified,
#' and vice-versa.
#' @return Just see an example, but here's a description:
#' Character matrix with entries starting with that in \code{copmat},
#' followed by the corresponding parameter values in \code{cparmat}
#' separated by commas, housed in parentheses.
#' @examples
#' G <- AtoG(CopulaModel::Dvinearray(4))
#'
#' ## Example 1
#' rv <- rvine(G, "frk", 2)
#' combine_copmat(rv)
#' combine_copmat(copmat=rv$copmat, cparmat=rv$cparmat)
#'
#' ## Example 2
#' rv <- rvine(G, "bvtcop", c(0.4, 5))
#' combine_copmat(rv)
#'
#' ## Example 3
#' rv <- rvine(matrix(4:1, ncol = 4))
#' combine_copmat(rv)
#' @export
combine_copmat <- function(rv, copmat, cparmat, digits=3) {
    if (!missing(rv)) {
        summat <- rv$copmat
        cparmat <- rv$cparmat
    } else {
        summat <- copmat
    }
    ## 1. Make parameter matrix a character matrix.
    cparmat <- apply(cparmat, 1:2, function(l) {
        l <- l[[1]]
        if (is.null(l)) {
            res <- ""
        } else {
            res <- paste(signif(l, digits), collapse=", ")
        }
        if (res != "") res <- paste0("(", res, ")")
        res
    })
    ## 2. Combine matrices
    for (i in seq_len(nrow(cparmat))) for (j in seq_len(ncol(cparmat))) {
        summat[i, j] <- paste0(summat[i, j], cparmat[i, j])
    }
    summat
}

#' @export
copmat.rvine <- function(rv) {
    rv$copmat
}

#' @export
copmat <- function(...) UseMethod("copmat")

#' @export
cparmat.rvine <- function(rv) {
    rv$cparmat
}

#' @export
cparmat <- function(...) UseMethod("cparmat")
