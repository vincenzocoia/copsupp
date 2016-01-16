#' Plot Copula Fits
#'
#' Plots normal scores plots, with fitted quantile curves.
#'
#' @param u Uniform data
#' @param v Uniform data
#' @param cops Vector of copula families to show fit for.
#' @param cpars List of copula parameters for the families in \code{cops}.
#' @param tauset Vector of quantile indices to plot.
#' @import ggplot2 CopulaModel plyr reshape2
#' @examples
#' set.seed(1234)
#' ## Get data
#' u <- runif(100)
#' v <- CopulaModel::qcondgum(runif(100), u, 4)
#'
#' ## Fit models:
#' cops <- c("frk", "gum", "joe", "mtcj")
#' cpars <- lapply(cops, function(cop_)
#' fitbicop_lh(u, v, families=cop_)$cpar)
#'
#' ## Get plots:
#' plotcop_qcurve(u, v, cops, cpars)
#'
#' ## How about on a validation set?
#' u <- runif(100)
#' v <- CopulaModel::qcondgum(runif(100), u, 4)
#' plotcop_qcurve(u, v, cops, cpars)
#'
#' ## Over the entire range of quantile functions?
#' ##  You can add ggplot layers too:
#' library(ggplot2)
#' plotcop_qcurve(u, v, cops, cpars, tauset = 1:19/20) +
#'   labs(title = "Normal Scores Plot", x = "x", y = "y")
#' @export
plotcop_qcurve <- function(u, v, cops, cpars, tauset = space_taus(10)) {
    ## Since I want different curves to be plotted on different panes,
    ##  I'll have to make the curves myself.
    ## 1. Set-up the plotting data frame.
    x <- qnorm(u)
    y <- qnorm(v)
    ngroups <- length(cops)
    n <- length(u)
    cops <- factor(cops, levels=unique(cops), ordered = TRUE)
    scatdat <- data.frame(copula = rep(cops, each = n),
                          qnorm_u = rep(x, ngroups),
                          qnorm_v = rep(y, ngroups))
    ## 2. Make (x,y) coordinates for each function.
    grid <- seq(from = min(x), to = max(x), length.out = 100)
    ugrid <- pnorm(grid)
    qconddfs <- list()
    ## Get data frames of predictions for each copula
    for (i in 1:length(cops)) {
        ## Quantile function of the copula:
        qcond <- get(paste0("qcond", cops[i]))
        ## Matrix of predictions over the grid:
        qcondfit <- sapply(tauset, function(tau_) {
            qnorm(qcond(tau_, ugrid, cpars[[i]]))
        })
        ## Melt the matrix
        colnames(qcondfit) <- tauset
        qcondfit <- melt(qcondfit)
        qcondfit$Var1 <- NULL  # This is just the row number of the original matrix.
        names(qcondfit) <- c("tau", "prediction")
        qcondfit$grid <- rep(grid, length(tauset))
        qcondfit$copula <- cops[i]
        qconddfs[[i]] <- qcondfit
    }
    ## Combine predictions from each copula into one data frame.
    qconddfs <- do.call(rbind, qconddfs)
    ## Output the plotting command
    ggplot(scatdat, aes(qnorm_u, qnorm_v)) +
        geom_point() +
        geom_line(aes(grid, prediction, colour = tau, group = tau), data = qconddfs) +
        facet_wrap(~ copula)
}

#' Plot CNQR scores
#'
#' Plots scores of different model in descending order for easy comparison.
#'
#' @param cops Vector of names of the models.
#' @param scores Vector of numeric scores of the models.
#' @return A plot (\code{ggplot} command) of the models' scores
#' in descending order. The labels are "copula" on the x-axis, and
#' "CNQR score" on the y-axis, but you can change that by literally adding
#' a layer \code{labs(x="Your label here", y="Your label here")} to the output.
#' @examples
#' set.seed(12)
#' plotcop_score(LETTERS[1:5], rexp(5))
#' @import ggplot2
#' @export
plotcop_score <- function(cops, scores, comparewith = NULL) {
    ord <- order(scores)
    cops <- factor(cops, levels = cops[ord], ordered=TRUE)
    plotdat <- data.frame(copula = cops[ord], score = scores[ord])
    if (is.null(comparewith)) {
        comparelayer <- NULL
        ylabel <- labs(y = "CNQR score")
    } else {
        comparelayer <- geom_hline(yintercept = comparewith, linetype = "dotted")
        ylabel <- labs(y = "CNQR score (compared with previous score)")
    }
    ggplot(plotdat, aes(x = copula)) +
        comparelayer +
        geom_point(aes(y = score)) +
        geom_linerange(aes(ymin = 0, ymax = score)) +
        ylabel
}
