#' Extract Variables from a Regular Vine
#'
#' Extract variables from a regular vine, possibly truncated, in the order
#' of the vine array.
#'
#' @param rv A regular vine object.
#' @return Vector of vine variables.
#' @examples
#' rv <- rvine(AtoG(CopulaModel::Dvinearray(5)))
#' rv <- relabel(rv, c(5, 2, 4, 3, 1))
#' vars(rv)
#'
#' rv <- trunc(rv, 2)
#' vars(rv)
#' @export
vars.rvine <- function(rv) {
    G <- rv$G
    if (nrow(G) == 0) return(integer(0))
    rv$G[1, ]
}

#' @export
vars <- function(...) UseMethod("vars")
