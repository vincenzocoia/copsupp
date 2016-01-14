#' Reorder a Vine, If Possible
#'
#' These functions reorder a vine according to a user-specified order.
#' \code{reordervarray} reorders a vine array only, whereas
#' \code{reorder.rvine} reorders an "rvine" object.
#'
#' @param rv Object of type "rvine"
#' @param G Vine array matrix
#' @param ord Vector consisting of the variables in \code{G} or \code{rv},
#' in the order you want them to be in.
#' @return Returns \code{NULL} if the vine cannot be re-ordered in the
#' requested manner. Otherwise, the re-ordered vine array (\code{reordervarray})
#' or "rvine" object (\code{reorder.rvine}) is returned.
#' @examples
#' ## reordervarray:
#' G <- AtoG(CopulaModel::Dvinearray(5))
#' reordervarray(G, 5:1)
#' reordervarray(truncvarray(G, c(0, 1, 2, 1, 3)), c(3, 2, 4, 1, 5))
#' reordervarray(G, c(5, 1, 4, 2, 3)) # NULL
#' reordervarray(G, c(5:1, 87)) # Error
#'
#'
#' ## reorder.rvine:
#' copmat <- makevinemat("a", c("b","c"), letters[4:6], letters[7:10], zerocol=T)
#' cparmat <- makevinemat(1, 2:3, 4:6, 7:10, zerocol=T)
#' rv <- rvine(G, copmat, cparmat)
#' reorder(rv, 5:1)
#' reorder(rv, c(3, 2, 4, 1, 5))
#' reorder(rv, c(5, 1, 4, 2, 3)) # NULL
#' reorder(rv, c(5:1, 87)) # Error
#'
#' ## Trivial cases:
#' reordervarray(matrix(0,1,0), integer(0))
#' reorder(rvine(matrix(5)), 5)
#' @rdname reorder
#' @export
reordervarray <- function(G, ord) {
    ## Method: Build "Gnew" layer-by-layer by querying G.
    ncolG <- ncol(G)
    nrowG <- nrow(G)
    if (ncolG == 0) {
        if (length(ord) != 0) {
            stop(paste0("Can't reorder variable(s) ",
                        paste(ord, collapse=", "),
                        " on a blank vine array."))
        } else {
            return(G)
        }
    }
    ovars <- G[1, ]  # Original variables
    if (!setequal(union(ovars, ord), ovars) | length(ord) != length(ovars))
        stop("Requested variables do not match the variables in G.")

    if (nrowG == 1) return(matrix(ord, nrow=1))
    if (ncolG == 2) return(makevinemat(ord[1], ord[2:1]))
    ntrunc <- trunclevel(G)
    ## Initiate Gnew. After the initial checks, it'll have at least 3 columns
    ##  and 2 rows.
    Gnew <- matrix(0, nrow=nrowG, ncol=ncolG)
    Gnew[1, ] <- ord
    ## Go column-by-column, looking for the edges.
    for (jnew in 2:ncolG) {
        ## MOVE TO NEXT VARIABLE IN THE ORDER
        ## The variable of this new layer:
        thisvar <- ord[jnew]
        ## Candidate variables to link with "thisvar" (must be already in array)
        candvars <- ord[1:(jnew-1)]
        ## Now see which of these candidate variables link up:
        for (inew in 2:nrowG) {
            ## MOVE TO THE NEXT ROW
            ## Get variables conditioned on:
            condsetindices <- 1+seq_len(inew-2)
            condset <- Gnew[condsetindices, jnew]
            ## What does 'thisvar' pair up with, given condset?
            pairedvar <- NULL
            for (j in 1:ncolG) {
                ## QUERY G FOR A MATCH
                ## Check that the condset matches
                thiscondset <- G[condsetindices, j]
                if (setequal(thiscondset, condset)) {
                    ## Yes, this column has the appropriate condset. Does it have the
                    ##  appropriate pair?
                    pair_ <- G[c(1, inew), j]
                    includes_thisvar <- thisvar %in% pair_
                    includes_candvar <- any(sapply(candvars, function(v_) v_ %in% pair_))
                    if (includes_thisvar & includes_candvar) {
                        ## Yes, there is a pair. Take it and stop querying G.
                        pairedvar <- setdiff(pair_, thisvar)
                        break
                    }
                }
            }
            ## Was a match for this row found?
            if (is.null(pairedvar)) {
                ## No match was found. Move on to the next layer of G.
                break
            } else {
                ## A match was found. Put it in G, and move to the next row.
                Gnew[inew, jnew] <- pairedvar
            }
        }
    }
    ## Now that all the columns in Gnew are filled in as much as possible,
    ##  check that all the edges from G ended up in Gnew.
    ntruncnew <- trunclevel(Gnew)
    if (sum(ntruncnew) != sum(ntrunc)) {
        ## Not all the edges of G ended up in Gnew, meaning that the vine
        ##  cannot be ordered in the requested way. Return NULL.
        return(NULL)
    } else {
        ## Success -- all edges in G ended up in Gnew. Output Gnew.
        return(Gnew)
    }
}

#' @rdname reorder
#' @export
reorder.rvine <- function(rv, ord) {
    ## Reorder vine array
    G <- rv$G
    if (ncol(G) == 0 & length(ord) == 0) return(rv)
    Gnew <- reordervarray(G, ord)
    if (is.null(Gnew)) return(NULL)
    ## Reorder copula and cparmat
    copmatnew <- reformcopmat(rv$copmat, Gnew, G)
    cparmatnew <- reformcopmat(rv$cparmat, Gnew, G)
    ## Output
    rvine(Gnew, copmatnew, cparmatnew)
}

#' @export
reorder <- function(...) UseMethod("reorder")
