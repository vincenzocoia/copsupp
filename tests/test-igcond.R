# ## Test the igcondinv function.
# library(testthat)
# p <- 0:10/10
#
# k <- 5
# eta <- 0:10
# (chk <- igcond(igcondinv(p, k, eta), k, eta))
# expect_equal(p, chk)
#
# k <- 1.11
# eta <- 0:10
# (chk <- igcond(igcondinv(p, k, eta), k, eta))
# expect_equal(p, chk)
#
# ## Start to run into problems for smaller k:
# k <- 1.10
# eta <- 0:10
# (chk <- igcond(igcondinv(p, k, eta), k, eta))
# # expect_equal(p, chk)
#
# ## But no problem if eta is zero, in which case the function is 1/t:
# k <- 1.10
# eta <- 0
# (chk <- igcond(igcondinv(p, k, eta), k, eta))
#
# ## Large k, small theta
# k <- 500
# eta <- 1
# (chk <- igcond(igcondinv(p, k, eta), k, eta))
# expect_equal(p, chk)
#
# ## Small k, Large theta
# k <- 1.5
# eta <- 500
# (chk <- igcond(igcondinv(p, k, eta), k, eta))
# expect_equal(p, chk)
#
# ## Both large
# k <- 500
# eta <- 500
# (chk <- igcond(igcondinv(p, k, eta), k, eta))
# expect_equal(p, chk)
