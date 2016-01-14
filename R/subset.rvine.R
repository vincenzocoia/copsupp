#' Subset a Regular Vine
#'
#' Find the vine connecting a subset of variables from a bigger vine, if
#' it exists. \code{subsetvarray} only finds the subset of a vine array.
#'
#' @param rv A regular vine object.
#' @param G A vine array matrix.
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
#' subsetvarray(G, c(2, 4, 3))
#'
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
#' @rdname subsetvine
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


#' @rdname subsetvine
#' @export
subsetvarray <- function(G, select, justcheck = FALSE) {
    ## Extract info
    v <- G[1, ]
    k <- length(select)
    ## The trivial cases:
    unknownvars <- setdiff(select, v)
    if (length(unknownvars) > 0) {
        warning(paste("Variables not in vine were selected to subset:",
                      paste(unknownvars, collapse = ", ")))
        if (justcheck) return(FALSE) else return(NULL)
    }
    if (k == 0) {
        if (justcheck) return(TRUE) else return(matrix(ncol=0, nrow=0))
    }
    if (k == 1) {
        if (justcheck) return(TRUE) else return(matrix(select))
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
    Glayers <- list(Gsub[1, 1])
    for (i in 2:k) Glayers[[i]] <- Gsub[1:nlink[i], i]
    do.call(makevinemat, Glayers)
}
