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
#' @import igraph CopulaModel VineCopula
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
    ## If there's no base vine, then make the base vine consist of the first
    ##  layer (...which itself should have no copulas).
    if (is.null(basevine)) {
        basevine <- rvine2fitrvine(dat, rvine(matrix(var[1])))
        var <- var[-1]
    }
    if (length(var) == 0) return(basevine)
    ## What variables are in the base vine?
    var_base <- vars(basevine)
    if (length(intersect(var, var_base)) != 0)
        stop(paste0("Variable(s) ", paste(intersect(var, var_base), collapse=", "),
                    " already exist in the base vine."))
    var_all <- c(var_base, var)
    if (max(var_all) > ncol(dat))
        stop(paste0("The variable(s) ",
                    paste(var_all[var_all>ncol(dat)], collapse=", "),
                    " refer to columns that do not exist in 'dat'."))
    ## We now have a fitted base vine with at least one variable, and at least one other
    ##  variable is being added. We first need to get an array, then choose
    ##  and fit copula models. We'll do so by adding one variable at a time.
    ## Setup:
    #### Truncation and vine array
    G_base <- basevine$G
    G <- layer$G
    if (is.null(G)) {
        ## Get truncation for each layer
        ntrunc <- layer$ntrunc
        maxtrunc <- length(var_base) + 1:length(var) - 1
        if (is.null(ntrunc)) ntrunc <- maxtrunc
        if (length(ntrunc) == 1) ntrunc <- rep(ntrunc, length(var))
        ntrunc <- pmin(ntrunc, maxtrunc)
        ## Setup array with only the nodes.
        G <- do.call(makevinemat, lapply(ntrunc, function(ntrunc_) rep(NA, ntrunc_+1)))
        # G[1, ] <- var
    } else {
        ntrunc <- trunclevel(G, overall = FALSE)
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


}



#' Fit a Vine to Data
#'
#' Fits a vine to data using \code{\link{gausstrvine.mst}} and then
#' fitting the vine array layer-by-layer using \code{\link{fitrvine_basic}}.
#' Only fits a regular vine, not a generalized vine.
#'
#' @param dat Data matrix with uniform scores.
#' @param vbls Vector of variables (column numbers in \code{dat}) to fit a model to.
#' @param ntrunc Truncation level of vine.
#' @param families Copula families to fit. Sorry, no 'indepcop'.
#' @examples
#' set.seed(123)
#' dat <- matrix(runif(500), ncol = 5)
#' fitrvine_basic(dat)
#' fitrvine_basic(dat, 3)
#' fitrvine_basic(dat, integer(0))
#' @import igraph CopulaModel VineCopula
#' @export
fitrvine_basic <- function(dat, vbls = 1:ncol(dat), ntrunc = length(vbls) - 1,
                           families = c("bvncop","bvtcop","mtcj","gum",
                                        "frk","joe","bb1","bb7","bb8")) {
    ## Temporary function -- don't allow for a basevine to be fit. Just fit everything.
    if (ntrunc == 0) {
        res <- rvine(matrix(vbls, nrow = 1))
        return(res)
    }
    k <- length(vbls)
    d <- k # Don't feel like changing it.
    if (k == 0) {
        res <- rvine(matrix(0,0,0))
        return(rvine2fitrvine(dat, res))
    }
    if (k == 1) {
        res <- rvine(matrix(vbls))
        return(rvine2fitrvine(dat, res))
    }
    if (k == 2) {
        v1 <- vbls[1]
        v2 <- vbls[2]
        res <- fitbicop_lh(dat[, v1], dat[, v2])
        res <- rvine(makevinemat(v1, c(v2, v1)),
                     copmat = res$cop, cparmat = res$cpar)
        return(rvine2fitrvine(dat, res))
    }
    ## Order data
    odat <- dat[, vbls]
    ## Get correlation matrix and choose vine array.
    cormat <- cor(qnorm(odat))
    library(igraph)  # Error if not present, even if igraph is under "import".
    A <- gausstrvine.mst(cormat, k-1)$RVM$VineA
#     G <- AtoG(A)
#     ## Fit the first two variables as the "base vine":
#     fit1 <- fitbicop_lh(odat[, A[1,1]], odat[, A[2,2]])
#     basevine <- rvine(makevinemat(A[1,1], A[2:1, 2]), fit1$cop, fit1$cpar)
#     ## Fit the next layers:
#     for (l in 3:k) {
#         e <- G[1:l, l]
#         thisfit <- fitlayer(odat, basevine, edges = e)
#         basevine <- addlayer(basevine, e, thisfit$cops, thisfit$cpars)
#     }
#     return(basevine)
    familyset <- sort(c(copname2num(families), recursive=T))
    names(familyset) <- NULL
    ## Now get and fit copulas. We'll use RVineCopSelect().
    ##  (use capture.output() to keep RVineCopSelect() quiet)
    capture.output(vinefit <- RVineCopSelect(odat,
                                             familyset = familyset,
                                             Matrix = A,
                                             trunclevel = ntrunc))
    ## Extract things.
    #### 1. copmat
    copmatind <- vinefit$family[(d:1)[1:ntrunc], d:1]
    if (!is.matrix(copmatind)) copmatind <- matrix(copmatind, ncol = d)
    copmat <- apply(copmatind, 1:2, copnum2name)
    copmat[!upper.tri(copmat)] <- ""
    #### cparmat
    parmat1 <- vinefit$par[(d:1)[1:ntrunc], d:1]
    parmat2 <- vinefit$par2[(d:1)[1:ntrunc], d:1]
    if (!is.matrix(parmat1)) parmat1 <- matrix(parmat1, ncol = d)
    if (!is.matrix(parmat2)) parmat2 <- matrix(parmat2, ncol = d)
    parvec <- numeric(0)
    len <- integer(0)
    for (i in 1:ntrunc) for (j in (i+1):d) {
        if (parmat1[i, j] == 0) {
            len <- c(len, 0)
        } else {
            parvec <- c(parvec, parmat1[i, j])
            if (parmat2[i, j] != 0) {
                parvec <- c(parvec, parmat2[i, j])
                len <- c(len, 2)
            } else {
                len <- c(len, 1)
            }
        }
    }
    cparmat <- makeuppertri.list(parvec, len, nrow = ntrunc, ncol = d)
    ## It seems that, when RVineCopSelect fits a 90- or 270-degree rotated
    ##  copula, it also makes the parameters negative. Need to fix this.
    ## joe; mtcj; gum; bb1; bb6; bb7; bb8
    for (i in 1:nrow(copmat)) for (j in (i+1):ncol(copmat)) {
        thiscop <- copmat[i, j]
        ## Is this copula model 90- or 270-degree rotated?
        lastchar <- substring(thiscop, first=nchar(thiscop), last=nchar(thiscop))
        if (lastchar == "u" | lastchar == "v") {
            ## Copula name ends in "u" or "v". Could that mean that it's "rotated"?
            wouldbe_pcop <- paste0("p", substring(thiscop, first=1, last=nchar(thiscop)-1))
            if (exists(wouldbe_pcop)) rotated <- TRUE else rotated <- FALSE
            ## Now that we know if it's rotated, flip the parameters if need be.
            if (rotated) cparmat[i, j][[1]] <- -cparmat[i, j][[1]]
        }
    }
    ## Return fit:
    G <- AtoG(A)[1:(ntrunc+1), ]
    newvbls <- G[1, ]
    G <- relabelvarray(G, vbls[newvbls])
    rvine(G, copmat, cparmat)
}

#' @export
print.fitrvine <- function(rv) {
    d <- ncol(rv$G)
    if (d == 0) return(cat("Empty fitted vine: no variables."))
    v <- vars(rv)
    ntrunc <- nrow(rv$G) - 1
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
