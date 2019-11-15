## Check that QYgX() gives the same result as rvineqcond()
library(CopulaModel)
library(copsupp)
library(cnqr)
library(testthat)
data(egdat)
xyinfo <- xylink(egvine)

## Useful quantities
A <- GtoA(egvine$G)
parvec <- c(t(egvine$cparmat), recursive=TRUE)
np <- apply(egvine$cparmat, 1:2, function(l) length(l[[1]]))
qcondmat <- apply(egvine$copmat, 1:2, function(cop) {
    if (cop != "") return(paste0("qcond", cop)) else return("")
})
pcondmat <- apply(egvine$copmat, 1:2, function(cop) {
    if (cop != "") return(paste0("pcond", cop)) else return("")
})

## Check pcondseq
uind_pcondseq <- pcondseq(egdat, ord=1:5, rv=egvine)
uind_rvinepcond <- rvinepcond(parvec, egdat, A=A, ntrunc=4, pcondmat=pcondmat, np=np)
expect_identical(uind_pcondseq, uind_rvinepcond)

## Compute quantiles
## NOTE: rvineqcond() computes the (5|1:4)-quantiles. That is, the quantile
##  levels used for the i'th observation is the i'th PIT score of 5|1:4.
## NOTE2: We need 'uind' here, separate from uind_pcondseq, because the
##  former has PIT scores (4), (3|4), (2|3:4), (1|2:4), which are needed in the
##  computation of the 5|1:4 quantiles. uind_pcondseq has PIT scores
##  (1), (2|1), (3|1:2), (4|1:3), which are not needed in the computation of
##  the 5|1:4 quantile.
uind <- pcondseq(egdat, ord=xyinfo$xord, rv=egvine)
yhat_QYgX <- sapply(1:length(vcond), function(i)
    QYgX(uind_pcondseq[i, 5], uind, cops=xyinfo$cops, cpars=xyinfo$cpars, QY=identity)[i, ])
yhat_rvineqcond <- rvineqcond(uind_pcondseq, A, ntrunc=4,
                              parvec=parvec,
                              np=np,
                              qcondmat=qcondmat,
                              pcondmat=pcondmat)[, 5]
expect_equal(yhat_QYgX, yhat_rvineqcond)  # They're equal, but not identical.
max(abs(yhat_QYgX - yhat_rvineqcond))
