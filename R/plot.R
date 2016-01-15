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
#' plotqcurve_edge(u, v, cops, cpars)
#'
#' ## How about on a validation set?
#' u <- runif(100)
#' v <- CopulaModel::qcondgum(runif(100), u, 4)
#' plotqcurve_edge(u, v, cops, cpars)
#'
#' ## Over the entire range of quantile functions?
#' ##  You can add ggplot layers too:
#' library(ggplot2)
#' plotqcurve_edge(u, v, cops, cpars, tauset = 1:19/20) +
#'   labs(title = "Normal Scores Plot", x = "x", y = "y")
#' @export
plotqcurve_edge <- function(u, v, cops, cpars, tauset = space_taus(10)) {
    ## Since I want different curves to be plotted on different panes,
    ##  I'll have to make the curves myself.
    ## 1. Set-up the plotting data frame.
    x <- qnorm(u)
    y <- qnorm(v)
    ngroups <- length(cops)
    n <- length(u)
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

# function(w = 1, g = identity)
