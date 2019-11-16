#' Incorporate Parameter Lengths
#'
#' Converts a vector of copula parameters, along with information about
#' how many copula parameters belongs to each family, into a list
#' of parameters for each family.
#'
#' @param cparvec Vector. Intended to be copula parameters.
#' @param len Vector of non-negative integers indicating how many parameters
#' belong to each family. Should have \code{sum(len)=length(cparvec)}.
#' @return List of copula parameters for each family corresponding to the
#' entries in \code{len}.
#' @examples
#' cparvec2cpar(-3:3/2, len = c(1, 2, 1, 0, 2, 1, 0))
#' @export
cparvec2cpar <- function(cparvec, len) {
    if (sum(len) != length(cparvec))
        stop("`len` indicates a different number of parameters than what's in `cparvec`.")
    cpar <- list()
    parnum <- 0
    for (i in 1:length(len)) {
        np <- len[i]
        cpar[[i]] <- cparvec[parnum + seq_len(np)]
        parnum <- parnum + np
    }
    return(cpar)
}
