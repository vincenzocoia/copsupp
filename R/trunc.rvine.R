#' Truncate a Regular Vine
#'
#' Truncates a vine (\code{trunc.rvine}), or a vine array (\code{truncvarray}).
#'
#' @param rv A regular vine object.
#' @param G A vine array.
#' @param ntrunc Integer; truncation level. Or, vector specifying the truncation
#' level of each column of the G-Vine array.
#' @return If \code{ntrunc >= nrow(G) - 1}, the original vine is
#' returned. Otherwise, an 'rvine' object is returned
#' with the truncation implemented.
#' @examples
#' rv <- rvine(AtoG(CopulaModel::Dvinearray(6)), "frk", 2)
#' trunc(rv, 3)
#' trunc(rv, 100)
#' trunc(rv, 0)
#'
#' ## Column-specific truncation
#' (rv1 <- trunc(rv, c(0, 1, 2, 1, 2, 4)))
#' lapply(rv1, identity)
#' rv2 <- trunc(rv, c(4, 4, 4, 1, 2, 4))
#' identical(rv1, rv2)
#' @rdname trunc
#' @export
trunc.rvine <- function(rv, ntrunc) {
    G <- truncvinemat(rv$G, ntrunc)
    copmat <- truncvinemat(rv$copmat, ntrunc, varray=FALSE, zero="")
    cparmat <- truncvinemat(rv$cparmat, ntrunc, varray=FALSE, zero=list(NULL))
    rvine(G, copmat, cparmat)
}

#' @export
trunc <- function(...) UseMethod("trunc")

#' @rdname trunc
#' @export
truncvarray <- function(G, ntrunc) {
    truncvinemat(G, ntrunc)
}

#' Get Truncation Level
#'
#' Extract the truncation level of a vine array. Intended for internal use.
#'
#' @param G Vine array.
#' @param overall Logical; \code{TRUE} returns the overall truncation level,
#' \code{FALSE} the truncation level of each column.
#' @examples
#' G <- AtoG(CopulaModel::Dvinearray(6))
#' G <- truncvinemat(G, c(0, 1, 2, 1, 2, 4))
#' trunclevel(G)
#' trunclevel(G, TRUE)
#' @export
trunclevel <- function(G, overall = FALSE) {
    r <- nrow(G)
    res <- r - apply(G, 2, function(col) sum(col == 0)) - 1
    if (overall) return(max(res)) else return(res)
}

#' Truncate a Vine Matrix
#'
#' Truncates a vine matrix (either the vine array, copmat, or cparmat).
#' Intended for internal use.
#'
#' @param mat Vine matrix
#' @param ntrunc Integer; truncation level. Or vector for each column.
#' @return A matrix, with the desired truncation.
#' @examples
#' (G <- AtoG(CopulaModel::Dvinearray(6)))
#' truncvinemat(G, 3)
#'
#' copmat <- makevinemat("frk", rep("frk", 2), rep("frk", 3), rep("frk", 4), zerocol=T)
#' truncvinemat(copmat, c(0, 1, 0, 3, 1), varray=FALSE, zero="")
#' @export
truncvinemat <- function(mat, ntrunc, varray = TRUE, zero = 0) {
    d <- ncol(mat)
    r <- nrow(mat)
    if (r == 0) return(mat)
    if (r == 1 & varray) return(mat)
    if (length(ntrunc) == 1) ntrunc <- rep(ntrunc, d)
    ntrunc <- pmin(ntrunc, 1:d-1)
    nzeroes <- r - ntrunc  # Number of "zeroes" to go after non-zeroes in each col
    if (varray) nzeroes <- nzeroes - 1
    for (j in 1:d) {
        if (varray) {
            mat[ntrunc[j] + 1 + seq_len(nzeroes[j]), j] <- zero
        } else {
            mat[ntrunc[j] + seq_len(nzeroes[j]), j] <- zero
        }
    }
    ## Remove unneccessary zero-rows:
    overalltrunc <- max(ntrunc)
    if (varray) {
        mat <- mat[1:(overalltrunc+1), ]
        if (!is.matrix(mat)) mat <- matrix(mat, nrow = 1)
    } else {
        mat <- mat[seq_len(overalltrunc), ]
        if (!is.matrix(mat)) {
            mat <- as.list(mat)[-1]
            mat <- do.call(makevinemat, c(mat, zerocol = TRUE))
        }
    }
    mat
}
