Notes written August 16, 2016.

Main functions are rvinepcond() and rvineqcond(),
for the conditional cdf and quantile function of the last variable
given the first (d-1).

rvineqcond() has inputs from rvinepcond() (applied to first d-1 variables)
so it is not efficient. But code is written to be modification of functions
in CopulaModel, and can be used to check against more efficient code.

Note that these functions are valid only if all bivariate copulas in the
vine are permutation symmetric. Some small changes are needed if there are 
permutation asymmetric bivariate copulas; see the top of hjcondvine.R. 

Hence the functions are not valid if gumu, gumv are used:
Gumbel copula for (U,V) -> (1-U,V) or (U,1-V) for negative dependence.

pcondvine.r: checks of rvinepcond() with direct calculations for
              5-dimensional D-vine and C-vine.
qcondvine.r: checks of rvineqcond() with direct calculations for
              5-dimensional D-vine and C-vine.
hjcondvine.R: rvinepcond() and rvineqcond(), file to be source("")
chkrvinecond.r: checks for functions in hjcondvine.R

