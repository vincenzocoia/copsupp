## Test the igcondinv function.
library(testthat)

p <- 0:10/10
k <- 5
eta <- 0:10
(chk <- igcond(igcondinv(p, k, eta), k, eta))
expect_equal(p, chk)

k <- 1.11
(chk <- igcond(igcondinv(p, k, eta), k, eta))
expect_equal(p, chk)

## Start to run into problems for smaller k:
k <- 1.10
(chk <- igcond(igcondinv(p, k, eta), k, eta))
# expect_equal(p, chk)

## But no problem if eta is zero, in which case the function is 1/t:
eta <- 0
(chk <- igcond(igcondinv(p, k, eta), k, eta))
