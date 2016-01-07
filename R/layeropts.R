#' Fitting options for a new Vine Layer
#'
#' When fitting a new layer(s) to a vine, use this function to specify
#' "known" components of the new layer(s), as well as
#'
#' @param ntrunc Truncation level. Could be a vector corresponding to the
#' truncation level for the variables \code{var} or \code{A[1, ]}.
#' @note The arrays here use the newer form, where variables go in row 1.
layeropts <- function(var=NULL, A=NULL, ntrunc=NULL, cops=NULL, cpar=NULL,
                   families = c("bvncop","bvtcop","mtcj","gum","frk","joe","bb1","bb7","bb8")){
    ## Deal with array-related things first.
    if (is.null(A) & is.null(var))
        stop("At least one of 'var' or 'A' must be specified.")
    if (!is.null(A) & !is.null(var)) {
        warning("Both 'var' and 'A' are specified -- ignoring 'var'.")
        var <- A[1, ]
    }
    if (is.null(A)) {
        ## Check that the var input is good.
        if (any(is.na(var)) | any(!is.numeric(var)))
            stop("'var' must be entirely integers, containing no NA's.")
        if (any(var%%1!=0)) # Check that they're integers. is.integer won't work.
            stop("'var' must be entirely integers, containing no NA's.")
    }
    list(var=var, A=A, ntrunc=ntrunc, copmat=cops, cparmat=cpar, families=families)
}
