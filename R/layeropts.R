#' Fitting options for a new Vine Layer
#'
#' When fitting a new layer(s) to a vine, use this function to specify
#' "known" components of the new layer(s), as well as
#'
#' @param ntrunc Truncation level. Could be a vector corresponding to the
#' truncation level for the variables \code{var} or \code{G[1, ]}.
#' @note The arrays here use the newer form, where variables go in row 1.
#' @return A list of partially-specified layers of a vine.
#'
#' Regarding the vine array info:
#'
#' \itemize{
#'      \item \code{$var} Vector of new variables that these layers add, not
#'      necessarily in order.
#'      \item \code{$ntrunc} Depending on how much information is input, could be
#'      \code{NULL}, an integer for maximum tree depth of these layers, or a
#'      vector of tree depth for each layer.
#'      \item \code{$G} Either \code{NULL}, or a vine array.
#' }
#'
#' Regarding copula and parameter info:
#'
#' \itemize{
#'      \item \code{$copmat} Copula matrix for these layers. No blank column
#'      to the left. Some entries may be \code{NA}.
#'      \item \code{$cparmat} Copula parameter matrix for these layers. No
#'      blank column to the left. Some entries may contain \code{NA}'s.
#' }
#' @export
layeropts <- function(var=NULL, G=NULL, ntrunc=NULL, cops=NULL, cpar=NULL,
                      families = c("bvncop","bvtcop","mtcj","gum","frk",
                                   "joe","bb1","bb7","bb8")){
    ## Deal with array-related things first.
    if (is.null(G) & is.null(var))
        stop("At least one of 'var' or 'G' must be specified.")
    if (!is.null(G) & !is.null(var)) {
        warning("Both 'var' and 'G' are specified -- ignoring 'var'.")
        var <- G[1, ]
    }
    if (is.null(G)) {  # In this case, var is specified, G is not.
        ## Check that the var input contains integers.
        if (any(is.na(var)) | any(!is.numeric(var)))
            stop("'var' must be entirely integers, containing no NA's.")
        if (any(var%%1!=0)) # Check that they're integers. is.integer won't work.
            stop("'var' must be entirely integers, containing no NA's.")
    }
    if (is.null(var)) {  # In this case, G is specified, var is not.
        var <- G[1, ]
    }
    list(var=var, G=G, ntrunc=ntrunc, copmat=cops, cparmat=cpar, families=families)
}
