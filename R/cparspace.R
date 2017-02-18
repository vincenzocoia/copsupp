#' Get Copula Parameter Bounds
#'
#' For a selection of bivariate copula families, returns the parameter space.
#'
#' @param cop Single copula family name (character).
#' @param fn Logical; should a function of the parameter space be returned? If
#' not, returns a list of the lower and upper bounds.
#' @return Returns a function that accepts a vector of copula
#' parameters (corresponding to \code{cop}) and returns \code{TRUE} if that
#' parameter is in the parameter space and \code{FALSE} otherwise.
#'
#' In the case that you enter a copula family that is not on record in this
#' function, a function
#' that always returns \code{TRUE} is returned.
#' @note
#' I didn't notice the \code{\link{cparbound}} function in the
#' \code{CopulaModel} package. A wrapper is still useful though, so as
#' to include reflected copulas (-u, -v, or -r), and to output a function
#' of the parameter space. Coming soon.
#' @examples
#' cparspace(c("frk", "bvtcop", "bvtcop", "gum"))
#' cparspace("bvncop")(0.8)
#' cparspace("thiscopuladoesnotexist")("89")
#' @export
cparspace <- function(cop, fn = TRUE) {
    ## Get default bounds:
    bnds <- list(frk = list(left = -30, right = 30),
                 gum = list(left = 1, right = 20),
                 gumu = list(left = 1, right = 20),
                 gumv = list(left = 1, right = 20),
                 gumr = list(left = 1, right = 20),
                 indepcop = list(left = numeric(0), right = numeric(0)),
                 bvtcop = list(left = c(-1, 0), right = c(1, Inf)),
                 bvncop = list(left = -1, right =1),
                 joe = list(left = 1, right = Inf),
                 joeu = list(left = 1, right = Inf),
                 joev = list(left = 1, right = Inf),
                 joer = list(left = 1, right = Inf),
                 mtcj = list(left = -1, right = Inf),
                 mtcju = list(left = -1, right = Inf),
                 mtcjv = list(left = -1, right = Inf),
                 mtcjr = list(left = -1, right = Inf),
                 bb1 = list(left = c(0, 1), right = c(Inf, Inf)),
                 bb1u = list(left = c(0, 1), right = c(Inf, Inf)),
                 bb1v = list(left = c(0, 1), right = c(Inf, Inf)),
                 bb1r = list(left = c(0, 1), right = c(Inf, Inf)),
                 bb6 = list(left = c(1, 1), right = c(Inf, Inf)),
                 bb6u = list(left = c(1, 1), right = c(Inf, Inf)),
                 bb6v = list(left = c(1, 1), right = c(Inf, Inf)),
                 bb6r = list(left = c(1, 1), right = c(Inf, Inf)),
                 bb7 = list(left = c(1, 0), right = c(Inf, Inf)),
                 bb7u = list(left = c(1, 0), right = c(Inf, Inf)),
                 bb7v = list(left = c(1, 0), right = c(Inf, Inf)),
                 bb7r = list(left = c(1, 0), right = c(Inf, Inf)),
                 bb8 = list(left = c(1, 0), right = c(Inf, 1)),
                 bb8u = list(left = c(1, 0), right = c(Inf, 1)),
                 bb8v = list(left = c(1, 0), right = c(Inf, 1)),
                 bb8r = list(left = c(1, 0), right = c(Inf, 1)),
                 bskewncop = list(left = -c(1,1,1), right=c(1,1,1)),
                 bskewncopp = list(left = -c(1,1,1), right=c(1,1,1)),
                 igcop = list(left = c(0, 1), right=c(Inf, Inf)),
                 igcopu = list(left = c(0, 1), right=c(Inf, Inf)),
                 igcopv = list(left = c(0, 1), right=c(Inf, Inf)),
                 igcopr = list(left = c(0, 1), right=c(Inf, Inf)),
                 iglcop = list(left=1, right=Inf),
                 iglcopu = list(left=1, right=Inf),
                 iglcopv = list(left=1, right=Inf),
                 iglcopr = list(left=1, right=Inf),
                 bb1rsk = list(left = c(0, 1, 0), right = c(Inf, Inf, 1)),
                 bb1sk = list(left = c(0, 1, 0), right = c(Inf, Inf, 1)),
                 bb1vsk = list(left = c(0, 1, 0), right = c(Inf, Inf, 1)),
                 bb1skp = list(left = c(0, 1, 0), right = c(Inf, Inf, 1)))
    ## Extract bounds:
    thesebounds <- bnds[cop]
    absent <- sapply(thesebounds, is.null)
    if (any(absent)) {
        mssng <- paste(unique(cop[absent]), collapse = ", ")
        warning(paste0("I don't have bounds for these copula families: ",
                       mssng, ". The parameter space function outputted will always",
                       " return TRUE."))
        return(function(cpar) TRUE)
    }
    left <- do.call(c, lapply(thesebounds, function(l) l$left))
    names(left) <- NULL
    right <- do.call(c, lapply(thesebounds, function(l) l$right))
    names(right) <- NULL
    if (fn) {
        if (cop %in% c("bskewncop", "bskewncopp"))
            return(function(cpar) all(cpar > left) & all(cpar < right) &
                       1-cpar[2]^2-cpar[3]^2-cpar[1]^2+2*cpar[2]*cpar[3]*cpar[1]>0)
        return(function(cpar) all(cpar > left) & all(cpar < right))
    } else {
        if (cop %in% c("bskewncop", "bskewncopp"))
            warning(paste("The parameter space for the skew normal copula is",
                          "not rectangular, and is a subset of this rectangle."))
        return(list(lower=left, upper=right))
    }
}
