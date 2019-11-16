# ## Comparison of copreg::rVineTruncCondCDF and copsupp::pcondD for D-Vines
# library(copreg)
# library(copsupp)
# library(CopulaModel)
#
# ## Make a vine model
# A <- truncvarray(Dvinearray(4), 2)
# copmat <- makeuppertri(c("gum", "mtcj", "bvncop",
#                          "frk", "new"), 2, 4, "")
# cparmat <- makeuppertri(c(2, 3, 0.8,
#                           2.5, 3.5), 2, 4)
# udat <- c(0.3, 0.4, 0.5, 0.6)
#
# ## Convert to other types of input:
# parvec <- t(cparmat)[lower.tri(t(cparmat))]
# pcondmat <- apply(copmat, 1:2, function(t) paste0("pcond", t))
#
# ## Evaluate
# rVineTruncCondCDF(parvec = parvec,
#                   A = A,
#                   udat = matrix(udat, nrow = 1),
#                   ntrunc = 2,
#                   pcondmat = pcondmat,
#                   np = makeuppertri(1, 2, 4))
# pcondD.generic(x = udat, i = 4, num = -3, copmat = copmat, cparmat = cparmat)
#
# ## Truncate A again:
# Anew <- truncvarray(A, 1)
# copmat <- reform.copmat(copmat, Anew, A)
# cparmat <- reform.copmat(cparmat, Anew, A)
# parvec <- t(cparmat)[lower.tri(t(cparmat))]
# pcondmat <- apply(copmat, 1:2, function(t) paste0("pcond", t))
# #### Inflate Anew -- I don't think rVineTruncCondCDF likes it when a vine array
# ####  is not either d-1 or d-2 truncated.
# Anew[2, ] <- 0
# Anew <- rbind(Anew, matrix(0, nrow = 2, ncol = 4))
# diag(Anew) <- 1:4
#
# ## Evaluate
# rVineTruncCondCDF(parvec = parvec,
#                   A = Anew,
#                   udat = matrix(udat, nrow = 1),
#                   ntrunc = 1,
#                   pcondmat = pcondmat,
#                   np = makeuppertri(1, 2, 4))
# pcondD.generic(x = udat, i = 4, num = -3, copmat = copmat, cparmat = cparmat)
#
# copmat <- rbind(copmat, cbind(matrix("", ncol = 1, nrow = 3),
#                               makeuppertri("indepcop", 3, 3, "")))
# cparmat <- makeuppertri.list(parvec, c(1, 1, 1, 0,0,0), 4, 4)
# pcondD.generic(x = udat, i = 4, num = -3, copmat = copmat, cparmat = cparmat)
