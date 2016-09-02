## Check that the conditional IG copula distributions work.
library(plyr)
library(copsupp)
library(testthat)
set.seed(123)

## Set up parameters
theta <- c(0.1, 1, 10, 100)
k <- c(1.1, 2.1, 5, 50, 100)
par <- expand.grid(theta=theta, k=k)

## Generate data
u <- runif(1000)
tau <- runif(1000)
d_ply(par, names(par), function(dfrow) {
    print(dfrow)
    theta_ <- dfrow$theta
    k_ <- dfrow$k
    v <- qcondigcop(tau, u, c(theta_, k_))
    tau2 <- pcondigcop(v, u, c(theta_, k_))
    # expect_that(tau, equals(tau2))
    print(plot(tau, tau2, main=paste(theta_, k_)))
})
