## Check that the igcop_helper_inv is working properly.
## If running from outside of the copsupp package, be sure to source()
##  the igcop_helper.R file.
library(CopulaModel)
library(testthat)

tau <- 0.5
chk1 <- tau

theta <- 10
k <- 5
(chk2 <- igcop_helper(igcop_helper_inv(tau, theta, k), theta, k))
expect_equal(chk1, chk2)

theta <- 10
k <- 10
(chk2 <- igcop_helper(igcop_helper_inv(tau, theta, k), theta, k))
expect_equal(chk1, chk2)

theta <- 0.4
k <- 10
(chk2 <- igcop_helper(igcop_helper_inv(tau, theta, k), theta, k))
expect_equal(chk1, chk2)
