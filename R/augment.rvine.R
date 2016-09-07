#' Add to a Column of an \code{rvine}
#'
#' Adds extra layers to a specified column of an \code{rvine} object,
#' including adding a new column altogether.
#'
#' @param obj The object of type \code{rvine} with which to augment.
#' @param a Vector of vine array indices to append to the vine array column.
#' @param cop Vector of copula families to append to the copula matrix column.
#' @param cpar List of copula parameters to append to the parameter matrix column.
#' @param col The column of the vine to add to (i.e., the common column
#' of the vine array, and copula and parameter matrices). \code{NULL} (default)
#' if you want to add a new column to the vine.
#' @note If you're adding a new column with only one variable, there's no
#' need to specify the \code{cop} and \code{cpar} arguments.
#'
#' Some checks are built-in to ensure that the resulting vine is indeed an
#' \code{rvine}, but it's not comprehensive. Would be nice to have some
#' function like \code{is.rvine} to check for sure.
#' @examples
#' ## Add to an empty vine:
#' rv <- rvine(matrix(nrow=0, ncol=0))
#' summary(rv)
#' rv2 <- augment(rv, a=4, col=1)
#' summary(rv2)
#'
#' ## Add to an independence vine:
#' rv <- rvine(matrix(4:1, ncol = 4))
#' summary(rv)
#' rv2 <- augment(rv, a=c(4,3),
#'                cop=c("bvtcop", "bvncop"),
#'                cpar=list(c(0.5, 3), -0.7),
#'                col=3)
#' summary(rv2)
#' summary(augment(rv2, integer(0)))
#'
#' ## You can't do some illegal things.
#' \dontrun{
#' augment(rv, a=1)
#' augment(rv, a=3:4, cop=c("frk", "frk"),
#'         cpar=list(4,3), col=2)
#' augment(rv, a=5:1, cop="frk", cpar=list(4))
#' }
#' ## But you can do some:
#' rv2 <- augment(rv, a=c(5,5), cop="frk", cpar=list(4))
#' summary(rv2)
#' @export
augment.rvine <- function(obj, a, cop, cpar, col=NULL) {
    if (length(a) == 0) return(obj)
    ## Extract the items
    G <- obj$G
    copmat <- obj$copmat
    cparmat <- obj$cparmat
    ## Append.
    ## First, check whether a new column is being requested to add.
    if (is.null(col)) col <- ncol(G) + 1
    if (col == ncol(G) + 1) {
        ## Yes. Add a column of blanks.
        ## Note: use `rep(0, nrow(G))` instead of just `0` because there might
        ##   be no rows in the first place.
        nrowG <- nrow(G)
        G       <- cbind(G,       rep(0,          nrowG))
        copmat  <- cbind(copmat,  rep("",         nrowG-1))
        cparmat <- cbind(cparmat, rep(list(NULL), nrowG-1))
        ## Just do a check: since a new column is being added, there should
        ##  be one less entry in `cop` and `cpar` than there is in `a`. If not,
        ##  throw an error because an invalid vine will result. But first,
        ##  since we want to allow for non-entry of cop and cpar in case 'a'
        ##  is of length 1, make these objects exist:
        if (missing(cop)) cop <- character(0)
        if (missing(cpar)) cpar <- list()
        if (length(cop) != length(a) - 1 | length(cpar) != length(a) - 1)
            stop(paste("Length of cop or cpar is not appropriate, and does not",
                       "result in a valid object of type 'rvine'."))
        ## Just another check: this new variable being added shouldn't
        ##  appear in the existing vine array.
        if (a[1] %in% G)
            stop(paste("Variable", a[1], "already exists in the vine, so it",
                       "can't start a new column."))
    } else {
        ## No, a new column is not to be added. This means a new variable is
        ##  not being introduced to the vine, and `a`, `cop`, and `cpar`
        ##  should all be of the same length. Just do a quick check.
        if (length(cop) != length(a) | length(cpar) != length(a))
            stop(paste("Length of cop or cpar is not appropriate, and does not",
                       "result in a valid object of type 'rvine'."))
    }
    ## Second, check whether there are enough rows in the matrices.
    #### How many zero (blank) entries are in the vine array column?
    zeroes <- sum(G[, col] == 0)
    nonzeroes <- nrow(G) - zeroes
    #### How many additional rows need to be added to accomodate the augmentation?
    addrows <- length(a) - zeroes
    #### Is this even "legal"? Can't exceed the matrix diagonal.
    if (nonzeroes + length(a) > col)
        stop(paste("Can't add that many rows to column", col, "of the rvine."))
    if (addrows > 0) {
        ## Not enough rows. Add them, with blanks.
        ncolG <- ncol(G)
        G       <- rbind(G,       matrix(0,          nrow=addrows, ncol=ncolG))
        copmat  <- rbind(copmat,  matrix("",         nrow=addrows, ncol=ncolG))
        cparmat <- rbind(cparmat, matrix(list(NULL), nrow=addrows, ncol=ncolG))
    }
    ## Finally, just go ahead and put in the new addition. There's room now.
    #### Indices in vine array to replace:
    G_ind <- min(which(G[, col] == 0)) + 1:length(a) - 1
    #### Indices in cparmat and copmat to replace. Should be the same as G_ind,
    ####  but shifted down by 1. This could result in the first index being 0
    ####  in the case where a new column is being appended, so just remove the
    ####  0 in that case.
    c_ind <- G_ind - 1
    if (c_ind[1] == 0) c_ind <- c_ind[-1]
    #### Do the replacement:
    G[G_ind, col] <- a
    if (length(c_ind) > 0) {
        copmat[c_ind, col] <- cop
        cparmat[c_ind, col] <- cpar
    }
    ## Modify the original object, and return it.
    ## Note: a new rvine object is not created from scratch, just in case
    ##   there are additional things included in the object that should be
    ##   preserved (such as with a 'cnqr' object in the `cnqr` package.)
    obj$G <- G
    obj$copmat <- copmat
    obj$cparmat <- cparmat
    return(obj)
}


#' @export
augment <- function(...) UseMethod("augment")
