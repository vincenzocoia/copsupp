#' Subset a Regular Vine
#'
#' Find the vine connecting a subset of variables from a bigger vine, if
#' it exists.
#'
#' @param rv A regular vine object.
#' @param select Vector of variables to subset,
#' if possible. The order of the variables does not matter.
#' @param justcheck Logical; should this function only check whether or not
#' the subset exists? \code{TRUE} if so.
#' @details Just a technicality:
#' by saying a subset "doesn't have an existing vine", I mean that
#' a vine can't be formed using nodes and edges from
#' the original -- not that the
#' joint distribution of the selected variables can't be created from a vine
#' (so as to say, for example, that the simplifying assumption of vines
#' doesn't hold for this distribution).
#' @return
#' If \code{justcheck} is \code{TRUE}, returns \code{TRUE} if the requested
#' subset exists, and \code{FALSE} if not.
#'
#' If \code{justcheck} is \code{FALSE}, returns
#' a vine of the subsetted variables, with variables
#' ordered according to their order in \code{G}; or \code{NULL} if
#' the subset does not form a vine.
#' @examples
#' ## Setup a vine.
#' G <- AtoG(CopulaModel::Dvinearray(5))
#' copmat <- makeuppertri(c("gum", "mtcj", "gal", "joe",
#'                          "frk", "gum", "bb7",
#'                          "bb1", "indepcop",
#'                          "bb8"), 4, 5, "")
#' cparmat <- makeuppertri.list(c(3, 2.5, 2, 1.5,
#'                                1, 1.3, 2, 2,
#'                                3, 4,
#'                                5, 0.5),
#'                                len = c(1,1,1,1,1,1,2,2,0,2),
#'                                4, 5)
#' (rv <- rvine(G, copmat, cparmat))
#'
#' ## Subset some variables.
#' subset(rv, c(2, 4, 3))
#' subset(rv, 5)
#' subset(rv, integer(0))
#'
#' ## This subset won't work:
#' subset(rv, c(4, 1), justcheck = TRUE)
#' ## But it will in a 0-truncated vine:
#' subset(trunc(rv, 0), c(4, 1), justcheck = TRUE)
#' subset(trunc(rv, 0), c(4, 1))
#'
#' ## Select variables not present?
#' subset(rv, c(2, 4, 17))
#' @export
subset.rvine <- function(rv, select, justcheck = FALSE) {
    ## Extract info
    G <- rv$G
    copmat <- rv$copmat
    cparmat <- rv$cparmat
    v <- vars(rv)
    k <- length(select)
    ## The trivial cases:
    unknownvars <- setdiff(select, v)
    if (length(unknownvars) > 0) {
        warning(paste("Variables not in vine were selected to subset:",
                   paste(unknownvars, collapse = ", ")))
        if (justcheck) return(FALSE) else return(NULL)
    }
    if (k == 0) {
        if (justcheck) return(TRUE) else return(rvine(matrix(ncol=0, nrow=0)))
    }
    if (k == 1) {
        if (justcheck) return(TRUE) else return(rvine(matrix(select)))
    }
    ## Get columns of G associated to the selection.
    colsel <- sort(sapply(select, function(s) which(v == s)))
    ## Take those columns of the vine array:
    Gsub <- G[, colsel]
    if (!is.matrix(Gsub)) Gsub <- matrix(Gsub, nrow = 1)
    nlink <- integer(0)
    ## Find out how much to subset for each column:
    for (i in 1:k) {
        ## For the i'th column, there can be anywhere from 1 to i consecutive
        ##  variables in 'select'.
        inselect <- Gsub[, i] %in% select
        thisnlink <- sum(inselect)
        ## Make sure the TRUEs are consecutive:
        if (sum(inselect[1:thisnlink]) < thisnlink) {
            if (justcheck) return(FALSE) else return(NULL)
        }
        nlink[i] <- thisnlink
    }
    if (justcheck) return(TRUE)
    ## Subset the matrices:
    copmatsub <- copmat[, colsel]
    if (nrow(copmat) == 1) copmatsub <- t(copmatsub)
    cparmatsub <- cparmat[, colsel]
    if (nrow(cparmat) == 1) cparmatsub <- t(cparmatsub)
    Glayers <- list(Gsub[1, 1])
    coplayers <- list()
    cparlayers <- list()
    for (i in 2:k) {
        Glayers[[i]] <- Gsub[1:nlink[i], i]
        coplayers[[i-1]] <- copmatsub[seq_len(nlink[i]-1), i]
        cparlayers[[i-1]] <- cparmatsub[seq_len(nlink[i]-1), i]
    }
    G <- do.call(makevinemat, Glayers)
    copmat <- do.call(makevinemat, c(coplayers, zerocol = TRUE))
    cparmat <- do.call(makevinemat, c(cparlayers, zerocol = TRUE))
    rvine(G, copmat, cparmat)
}

#' @export
subset <- function(...) UseMethod("subset")


#' Subset a vine array
#'
#' Find the vine connecting a subset of variables from a bigger vine, if
#' it exists. Deprecated; use \code{\link{subset.rvine}} instead.
#'
#' @param A vine array
#' @param select Vector of variable indices in \code{diag(A)} to subset,
#' if possible. The order of the variables does not matter.
#' @details Just a technicality:
#' by saying a subset "doesn't have an existing vine", I mean that
#' a vine can't be formed using nodes and edges from
#' the original -- not that the
#' joint distribution of the selected variables can't be created from a vine
#' (so as to say, for example, that the simplifying assumption of vines
#' doesn't hold for this distribution).
#' @return Returns a vine array of the subsetted variables, with variables
#' ordered according to their order in \code{A}; or \code{NULL} if
#' the subset does not form a vine.
#' @examples
#' ## Try a D-Vine
#' rv <- rvine(CopulaModel::Dvinearray(5))
#' subset(rv, c(2, 3, 4))
#' subset(rv, c(4, 1))
#' @export
rvinesubset <- function(A, select) {
    warning("'rvinesubset' is deprecated. Please use 'subset.rvine' instead.")
    p <- ncol(A)
    ntrunc <- nrow(A) - 1
    k <- length(select)
    if (k == p) return(A)
    if (k == 1) return(matrix(select))
    if (k == 0) return(A[integer(0), integer(0)])
    vars <- varray.vars(A)
    notselect <- setdiff(vars, select)
    ## Indices to keep (i.e. what columns of A, or which of the
    ##   ordered variables to keep?)
    ikeep <- sort(sapply(select, function(s) which(vars == s)))
    ilose <- setdiff(1:p, ikeep)
    ## Convert vine array to a "convenient" vine array, where the labels go on
    ##  top instead of the "skewed diagonal":
    Acon <- A[1:ntrunc, ]
    if (!is.matrix(Acon)) Acon <- matrix(Acon, nrow = ntrunc)
    diag(Acon) <- 0
    Acon <- rbind(matrix(vars, nrow = 1), Acon)
    ## Remove bottom rows if need be:
    if (k - 1 < ntrunc) Acon <- Acon[1:k, ]
    nrownew <- nrow(Acon)
    ## Remove columns of non-selected variables:
    Acon <- Acon[, ikeep]
    ## Go column-by-column, and gather variables that are in the selection set.
    unstrung <- integer(0)
    for (j in 1:ncol(Acon)) for (i in 1:nrow(Acon)) {
        if (Acon[i, j] %in% select) unstrung <- c(unstrung, Acon[i, j])
    }
    if (length(unstrung) != nrownew * k - choose(nrownew, 2)) return(NULL)
    res <- makeuppertri(unstrung, row = nrownew, col = k,
                        incDiag = TRUE, byRow = FALSE)
    ## Convert back to standard array format:
    varsnew <- res[1, ]
    res <- res[c(2:nrownew, 1), ]
    diag(res) <- varsnew[1:nrownew]
    res[nrownew, 1:(nrownew-1)] <- 0
    res
}
