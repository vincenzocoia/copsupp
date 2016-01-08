#' Fit a Vine Model
#'
#' Fits a joint distribution for the data using a vine copula. The vine array
#' is first chosen using the minimum spanning tree algorithm using
#' the function \code{gausstrvine.mst} in the \code{CopulaModel}
#' package, then the pairwise
#' copula models are chosen and fit using \code{RVineCopSelect} in the
#' \code{VineCopula} package.
#'
#' @param dat Matrix of data having Uniform marginal distributions;
#' columns represent variables, and rows observations.
#' @param rv Object of class "rvine", representing a pre-specification of the
#' vine to fit. Or, \code{NULL} for fully unspecified. See details.
#' @param var Vector of integers specifying the column numbers of \code{dat}
#' to fit a model to. Default is all variables.
#' @param ntrunc Integer, either \code{1, 2, ...,ncol(dat)-1},
#' of the truncation level of the vine to be fit.
#' @param families A vector of copula family names to try
#' fitting (will also consider their rotations/reflections). Limited to
#' those families available in \code{VineCopula} package, listed in
#' \code{\link{BiCopSelect}}.
#' @param ... Other arguments to pass to \code{VineCopula::RVineCopSelect}.
#' @return A "fitrvine" object, which has class \code{c("fitrvine", "rvine")},
#' which is a named list of the following:
#'
#' \itemize{
#'      \item \code{$A}: Vine array, truncated to \code{ntrunc}.
#'      \item \code{$copmat}: \code{ntrunc x ncol(A)} upper-triangular
#'      matrix of copula model names.
#'      \item \code{$cparmat}: \code{ntrunc x ncol(A)} upper-triangular
#'      matrix of copula parameters. Each entry is a list of length one containing
#'      the vector of copula parameters for that copula family.
#'      \item \code{$dat}: The inputted data matrix, \code{dat}.
#'      \item \code{$aic}: The AIC of the fitted model.
#'      \item \code{$bic}: The BIC of the fitted model.
#'      \item \code{$nllh}: The negative log likelihood of the fitted model.
#'      \item \code{$covmat}: Covariance matrix of the fitted parameters (in
#'      reading-order of the parameters, i.e.
#'      \code{c(t(cparmat)[lower.tri(t(cparmat))], recursive=TRUE)}).
#' }
#' @details
#' If you want to specify parts of the vine, then specify them in an "rvine" object
#' (see \code{\link{rvine}}):
#'
#' \enumerate{
#'      \item Your first option is to specify the vine array \code{A}.
#'      \item If the vine array is specified, then you can specify some or all
#'      of the copula families by putting them in \code{copmat}. Leave
#'      unspecified edges as \code{NA}.
#'      \item If there are copula families specified, you can specify parameters
#'      for those families by putting the parameters in \code{cparmat}.
#'      Unspecified parameters should be \code{NA}.
#' }
#'
#' For parts of the vine that are unspecified, you have some fitting options:
#'
#' \itemize{
#'      \item If you didn't specify a vine array \code{A}, you can select which
#'      variables (column numbers of \code{dat}) you'd like to fit through the
#'      \code{var} argument. You can also select the truncation level of
#'      the vine array through the argument \code{ntrunc}.
#'      \item If you left some copula families unspecified, you can indicate
#'      the candidate families in the \code{families} argument.
#' }
#' @import igraph, CopulaModel, VineCopula
#' @examples
#' ## Get some simulated data:
#' set.seed(152)
#' ntrunc <- 2
#' d <- 4
#' A0 <- truncvarray(CopulaModel::Dvinearray(d), ntrunc)
#' copmat0 <- makeuppertri("frk", ntrunc, d, "")
#' cparmat0 <- makeuppertri(3, ntrunc, d)
#' dat <- fvinesim(100, A0, copmat0, cparmat0)
#'
#' ## Fit a model to the data:
#' fit.rvine(dat, ntrunc=ntrunc)
#' fit.rvine(dat, c(4, 2, 3))
#' @export


fitrvine <- function(dat, layer=layeropts(1:ncol(dat)), basevine = NULL, ...) {
    ## Look at input
    if (is.data.frame(dat)) dat <- as.matrix(dat)
    n <- nrow(dat)
    var <- layer$var
    ## If no layers are being added, output the base vine.
    if (length(var) == 0) if (is.null(basevine)) {
        res <- rvine(matrix(nrow=0, ncol=0))
        return(rvine2fitrvine(dat, res))
    } else {
        if (is.fitrvine(basevine)) {
            return(basevine)
        } else {
            return(rvine2fitrvine(dat, basevine))
        }
    }
    ## If there's no base vine, then make the base vine consist of the first variable.
    if (is.null(basevine)) {
        basevine <- rvine2fitrvine(dat, rvine(matrix(var[1])))
        var <- var[-1]
    }
    if (length(var) == 0) return(basevine)
    ## What variables are in the base vine?
    var_base <- vars(basevine)
    if (length(intersect(var, var_base)) == 0)
        stop(paste0("Variables ", paste(intersect(var, var_base), collapse=", "),
                    " already exist in the base vine."))
    var_all <- c(var_base, var)
    if (max(var_all) > ncol(dat))
        stop(paste0("The variable(s) ",
                    paste(var_all[var_all>ncol(dat)], collapse=", "),
                    " refer to columns that do not exist in 'dat'."))
    ## We now have a fitted base vine with at least one variable, and at least one other
    ##  variable is being added. We first need to get an array, then chooose
    ##  and fit copula models. We'll do so by adding one variable at a time.
    ## Setup:
    #### Truncation and vine array
    A_base <- basevine$A
    A <- layer$A
    if (is.null(A)) {
        ntrunc <- layer$ntrunc
        maxtrunc <- length(var_base) + 1:length(var) - 1
        if (is.null(ntrunc)) ntrunc <- maxtrunc
        if (length(ntrunc) == 1) ntrunc <- rep(ntrunc, length(var))
        ntrunc <- pmin(ntrunc, maxtrunc)
        A <- do.call(makevinemat, lapply(ntrunc, function(ntrunc_) rep(NA, ntrunc_+1)))
        A[1, ] <- var
    } else {
        ntrunc <- apply(A, 2, function(col) nrow(A) - 1 - sum(col == 0, na.rm = TRUE))
    }
    #### Copula Matrix
    copmat <- layer$cops
    if (is.null(copmat)) {
        cop_layers <- lapply(ntrunc, function(ntrunc_) rep(list(layer$families), ntrunc_))
        copmat <- do.call(makevinemat, c(cop_layers, zerocol = TRUE))
    }
    copmatnew <- copmat
    #### Copula parameter matrix
    cparmat <- layer$cpars
    if (is.null(cparmat)) {
        cpar_layers <- lapply(ntrunc, function(ntrunc_) rep(list(NA), ntrunc_))
        cparmat <- do.call(makevinemat, c(cpar_layers, zerocol = TRUE))
    }
    cparmatnew <- cparmat
    ## Now add the layers.
    for (j in 1:length(var)) {
        var_ <- var[j]
        ## --- 1. Vine array ---
        ## In this section, we'll fill-in the missing parts
        ##  of the original vine array in the order of
        ##  highest partial correlation (approximated by lm()).
        A[1:(ntrunc[j]+1), j] <- fitseq(dat, A[1, j],
                                        linkwith = A[1, 1:j],
                                        a = A[1+1:ntrunc[j], j])
        ## We now have a completed column in A.
        ## --- 2. Fit the new layer ---
        ## Copula families
        fams <- copmat[1:ntrunc[j], j]
        if (!is.list(fams)) fams <- as.list(fams)
        fams <- lapply(fams, function(vec) if (all(is.na(vec))) layer$families else vec)
        ## Parameters
        speccpar <- cparmat[1:ntrunc[j], j]
        if (!is.list(speccpar)) speccpar <- as.list(speccpar)
        fitlayer_ <- fitlayer(dat, basevine, edges = A[1+0:ntrunc[j], j],
                              cops = fams, cpars = speccpar)
        ## --- 3. Update ---
        ## Update the base vine with the new layer
        cparmat[]
    }
    ## Re-estimate the layers altogether (they were originally done edge-by-edge)






#     ## The trivial case
#     if (d == 1 | d == 0) {
#         res <- rvine(matrix(var))
#         res <- c(res, dat=dat, aic=NA, bic=NA, nllh=NA, covmat=NA)
#         class(res) <- c("fitrvine", "rvine")
#         return(res)
#     }
#     if (ntrunc == 0) {
#         res <- rvine(matrix(var, nrow=1))
#         res <- c(res, dat=dat, aic=NA, bic=NA, nllh=NA, covmat=NA)
#         class(res) <- c("fitrvine", "rvine")
#         return(res)
#     }
#
#     ## Get specified things
#     if (is.null(rv)) {
#         ## Get vine array
#         if (d == 2) {
#             A <- matrix(c(1,0,1,2), ncol = 2) # Needs to have labels = 1:2.
#         } else {
#             ## Get correlation matrix
#             cormat <- cor(qnorm(dat[, var]))
#             ## Choose vine array
#             arrayfit <- gausstrvine.mst(cormat, ntrunc)
#             A <- arrayfit$RVM$VineA
#         }
#         ## Indicate blank copmat and cparmat
#         copmat <- NULL
#         cparmat <- NULL
#     } else {
#         ## Change labels to 1:d
#         A <- relabel(rv, 1:d)$A
#         ## Fill-in truncated part:
#         A <- rbind(A, matrix(0, nrow = d - ntrunc - 1, ncol = d))
#         diag(A) <- 1:d
#     }
#     familyset = sort(unique(c(copname2num(families), recursive = TRUE)))
#     ## Now get and fit copulas. We'll use RVineCopSelect() for now, which
#     ##  means we'll fit the entire vine first and then replace the fit
#     ##  with the pre-specified things. It's not ideal but it's a start.
#     ##  (use capture.output() to keep RVineCopSelect() quiet)
#     capture.output(vinefit <- RVineCopSelect(dat,
#                                              familyset = familyset,
#                                              Matrix = A,
#                                              trunclevel = ntrunc,
#                                              ...))
#     ## Extract things.
#     #### 1. copmat
#     copmatind <- vinefit$family[(d:1)[1:ntrunc], d:1]
#     if (!is.matrix(copmatind)) copmatind <- matrix(copmatind, ncol = d)
#     copmat <- apply(copmatind, 1:2, copnum2name)
#     copmat[!upper.tri(copmat)] <- ""
#     #### Replace copulas with specified copulas:
#     #### cparmat
#     parmat1 <- vinefit$par[(d:1)[1:ntrunc], d:1]
#     parmat2 <- vinefit$par2[(d:1)[1:ntrunc], d:1]
#     if (!is.matrix(parmat1)) parmat1 <- matrix(parmat1, ncol = d)
#     if (!is.matrix(parmat2)) parmat2 <- matrix(parmat2, ncol = d)
#     parvec <- numeric(0)
#     len <- integer(0)
#     for (i in 1:ntrunc) for (j in (i+1):d) {
#         if (parmat1[i, j] == 0) {
#             len <- c(len, 0)
#         } else {
#             parvec <- c(parvec, parmat1[i, j])
#             if (parmat2[i, j] != 0) {
#                 parvec <- c(parvec, parmat2[i, j])
#                 len <- c(len, 2)
#             } else {
#                 len <- c(len, 1)
#             }
#         }
#     }
#     if (all(len == 1)) {
#         cparmat <- makeuppertri(parvec, nrow = ntrunc, ncol = d)
#     } else {
#         cparmat <- makeuppertri.list(parvec, len, nrow = ntrunc, ncol = d)
#     }
#     ## It seems that, when RVineCopSelect fits a 90- or 270-degree rotated
#     ##  copula, it also makes the parameters negative. Need to fix this.
#     ## joe; mtcj; gum; bb1; bb6; bb7; bb8
#     for (i in 1:nrow(copmat)) for (j in (i+1):ncol(copmat)) {
#         thiscop <- copmat[i, j]
#         ## Is this copula model 90- or 270-degree rotated?
#         lastchar <- substring(thiscop, first=nchar(thiscop), last=nchar(thiscop))
#         if (lastchar == "u" | lastchar == "v") {
#             ## Copula name ends in "u" or "v". Could that mean that it's "rotated"?
#             wouldbe_pcop <- paste0("p", substring(thiscop, first=1, last=nchar(thiscop)-1))
#             if (exists(wouldbe_pcop)) rotated <- TRUE else rotated <- FALSE
#             ## Now that we know if it's rotated, flip the parameters if need be.
#             if (rotated) {
#                 if (is.list(cparmat[1,1])) {
#                     cparmat[i, j][[1]] <- -cparmat[i, j][[1]]
#                 } else {
#                     cparmat[i, j] <- -cparmat[i, j]
#                 }
#             }
#         }
#     }
#     ## Output results
#     Avars <- varray.vars(A)
#     list(cdf=margs,
#          A=relabel.varray(A, var[Avars]),
#          copmat=copmat,
#          cparmat=cparmat)
}

#' @export
print.fitrvine <- function(rv) {
    d <- ncol(rv$A)
    if (d == 0) return(cat("Empty fitted vine: no variables."))
    v <- var(rv)
    ntrunc <- nrow(rv$A) - 1
    if (ntrunc == 0) {
        trunctext <- "Independent 'fitted'"
    } else {
        trunctext <- paste0(ntrunc, "-truncated fitted")
    }
    cat(paste0(trunctext, " vine with variables ", paste(v, collapse = ", "), ".\n"))
    cat(paste0("aic,bic: ", rv$aic, ",", rv$bic, "\n"))
    invisible()
}

#' @export
is.fitrvine <- function(rv) {
    inherits(rv, "fitrvine")
}
