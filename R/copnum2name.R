#' Convert between Copula Number and Name
#'
#' The \code{VineCopula} package codes the copula families by some integer, whereas
#' the CopulaModel package uses a naming convention. These functions converts
#' the integer code to the copula family name (\code{copnum2name}), or
#' vice versa (\code{copname2num}).
#'
#' @param num Vector of integer codes.
#' @param name Vector of copula family names (without indication of rotation)
#' @return 
#' For \code{copnum2name}, returns a character vector of the copula family
#' names (down to the rotation). \code{NA} for non-existent families.
#' 
#' For \code{copname2num}, returns a list with entries corresponding to
#' \code{name} of the integer codes for that family, including all the rotations
#' if they're different. \code{NULL} entry if copula model not found.
#' @seealso See the help page for \code{\link{RVineCopSelect}}
#' in the \code{VineCopula} package for the mapping from code to copula family.
#' @rdname copnum2name
#' @export
copnum2name <- function(num) {
    x <- num + 1
    fams <- c("indepcop",
              "bvncop",
              "bvtcop",
              "mtcj",
              "gum",
              "frk",
              "joe",
              "bb1",
              "bb6",
              "bb7",
              "bb8", NA, NA,
              "mtcjr",
              "gumr", NA,
              "joer",
              "bb1r",
              "bb6r",
              "bb7r",
              "bb8r", NA, NA,
              "mtcju",
              "gumu", NA,
              "joeu",
              "bb1u",
              "bb6u",
              "bb7u",
              "bb8u", NA, NA,
              "mtcjv",
              "gumv", NA,
              "joev",
              "bb1v",
              "bb6v",
              "bb7v",
              "bb8v")
    fams[x]
}

#' @rdname copnum2name
#' @export
copname2num <- function(name) {
    f <- list(indepcop = 0,
              bvncop = 1,
              bvtcop = 2,
              mtcj = c(3, 13, 23, 33),
              gum = c(4, 14, 24, 34),
              frk = 5,
              joe = c(6, 16, 26, 36),
              bb1 = c(7, 17, 27, 37),
              bb6 = c(8, 18, 28, 38),
              bb7 = c(9, 19, 29, 39),
              bb8 = c(10, 20, 30, 40))
    res <- f[name]
    names(res) <- name  # Ensures NULL output is named appropriately.
    nulls <- sapply(res, is.null)
    if (any(nulls)) warning(paste0("Copula models could not be found in VineCopula pkg: ",
                                   paste(unique(name[nulls]), collapse = ", ")))
    res
}