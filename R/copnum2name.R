#' Convert Copula Number to Name
#'
#' The VineCopula package codes the copula families by some integer, whereas
#' the CopulaModel package uses a naming convention. This function converts
#' the integer code to the copula family name.
#'
#' @param num Integer code.
#' @return A string of the name of the copula family.
#' @seealso See the help page for \code{\link{RVineCopSelect}}
#' in the \code{VineCopula} package for the mapping from code to copula family.
#' @export
copnum2name <- function(num) {
    x <- num
    if (x == 0) return("indepcop")
    fams <- c("bvncop",
              "bvtcop",
              "mtcj",
              "gum",
              "frk",
              "joe",
              "bb1",
              "bb6",
              "bb7",
              "bb8", "", "",
              "mtcjr",
              "gumr", "",
              "joer",
              "bb1r",
              "bb6r",
              "bb7r",
              "bb8r", "", "",
              "mtcju",
              "gumu", "",
              "joeu",
              "bb1u",
              "bb6u",
              "bb7u",
              "bb8u", "", "",
              "mtcjv",
              "gumv", "",
              "joev",
              "bb1v",
              "bb6v",
              "bb7v",
              "bb8v")
    if (x > 40)
        stop("Are you asking for a Tawn copula family? There's no functionality for it yet.")
    fams[x]
}
