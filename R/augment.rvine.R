augment.rvine <- function(obj, col, a, cop, cpar) {
    if (length(a) == 0) return(obj)
    ## Extract the items
    G <- obj$G
    copmat <- obj$copmat
    cparmat <- obj$cparmat
    ## Append.
    ## First, check whether a new column is being requested to add.
    if (col == ncol(G) + 1) {
        ## Yes. Add a column of blanks.
        ## Note: use `rep(0, nrow(G))` instead of just `0` because there might
        ##   be no rows in the first place.
        nrowG <- nrow(G)
        G       <- cbind(G,       rep(0,          nrowG))
        copmat  <- cbind(copmat,  rep("",         nrowG))
        cparmat <- cbind(cparmat, rep(list(NULL), nrowG))
        ## Just do a check: since a new column is being added, there should
        ##  be one less entry in `cop` and `cpar` than there is in `a`. If not,
        ##  throw an error because an invalid vine will result.
        if (length(cop) != length(a) - 1 | length(cpar) != length(a) - 1)
            stop(paste("Requested augmentation in augment.rvine() does not",
                       "result in a valid object of type 'rvine'."))
    } else {
        ## No, a new column is not to be added. This means a new variable is
        ##  not being introduced to the vine, and `a`, `cop`, and `cpar`
        ##  should all be of the same length. Just do a quick check.
        if (length(cop) != length(a) | length(cpar) != length(a))
            stop(paste("Requested augmentation in augment.rvine() does not",
                       "result in a valid object of type 'rvine'."))
    }
    ## Second, check whether there are enough rows in the matrices.
    #### How many zero (blank) entries are in the vine array column?
    zeroes <- sum(G[, col] == 0)
    #### How many additional rows need to be added to accomodate the augmentation?
    addrows <- length(a) - zeroes
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
    copmat[c_ind, col] <- cop
    cparmat[c_ind, col] <- cpar
    ## Modify the original object, and return it.
    ## Note: a new rvine object is not created from scratch, just in case
    ##   there are additional things included in the object that should be
    ##   preserved (such as with a 'cnqr' object in the `cnqr` package.)
    obj$G <- G
    obj$copmat <- copmat
    obj$cparmat <- cparmat
    return(obj)
}
