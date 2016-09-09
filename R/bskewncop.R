## Bivariate skew normal copula: functions from Harry
## NOTE: I'm not going to add any more documentation than what's already here.

#' Skew Normal Copula Functions
#'
#' skew-normal copula from Azzalini-DallaValle's bivariate skew-normal
#' u,v = values or vectors in (0,1)
#' param = (rho,ga1,ga2), with -1<rho<1, -1<ga1<1, -1<ga2<1
#'  and a constraint 1-ga1^2-ga2^2-rho^2+2*ga1*ga2*rho>0
#'
#' param=c(rho,ga1,ga2), constraint is that
#' correlation matrix with these parameters should be positive definite
#' @rdname bskewncop
#' @import sn CopulaModel
#' @export
dbskewn=function(y1,y2,param)
{ rho=param[1]; ga1=param[2]; ga2=param[3]
r2=1-rho^2
s=1-ga1^2-ga2^2-rho^2+2*ga1*ga2*rho # should be positive
s=sqrt(s/r2)
be1=(ga1-rho*ga2)/r2
be2=(ga2-rho*ga1)/r2
bpdf=2*dbvn2(y1,y2,rho)*pnorm((be1*y1+be2*y2)/s)
bpdf
}

#============================================================


#' @return
#' \code{dbskewncop}: bivariate skew-normal copula density
#' @rdname bskewncop
#' @import sn CopulaModel
#' @export
dbskewncop=function(u,v,param)
{ rho=param[1]; ga1=param[2]; ga2=param[3]
r2=1-rho^2
a1=ga1/sqrt(1-ga1^2)
a2=ga2/sqrt(1-ga2^2)
s=1-ga1^2-ga2^2-rho^2+2*ga1*ga2*rho # should be positive
s=sqrt(s/r2)
be1=(ga1-rho*ga2)/r2
be2=(ga2-rho*ga1)/r2
y1=qsn(u,0,1,alpha=a1)
y2=qsn(v,0,1,alpha=a2)
bpdf=2*dbvn2(y1,y2,rho)*pnorm((be1*y1+be2*y2)/s)
denom=dsn(y1,alpha=a1)*dsn(y2,alpha=a2)
cpdf=bpdf/denom
cpdf
}

#' @return
#' \code{pbskewncop}: bivariate skew-normal copula cdf
#' @rdname bskewncop
#' @import sn CopulaModel
#' @export
pbskewncop=function(u,v,param)
{ rho = param[1]; gav = param[2:3]
rmat=matrix(c(1,rho,rho,1),2,2)
tem=solve(rmat,gav)
alpv=tem/sqrt(1-sum(gav*tem))
ga1=gav[1]; ga2=gav[2]
a1=ga1/sqrt(1-ga1^2)
a2=ga2/sqrt(1-ga2^2)
y1=qsn(u,alpha=a1); y2=qsn(v,alpha=a2)
cdf=pmsn(c(y1,y2),Omega=rmat,alpha=alpv)
cdf
}

#' @note  pcond21 = P(V<v|U=u) for bskewn copula with  y1=qsn(u), y2=qsn(v)
#' @return
#' \code{pcondbskewncop21}: conditional cdf P(V<v|U=u) of bivariate skew-normal copula
#' @rdname bskewncop
#' @import sn CopulaModel
#' @export
pcondbskewncop21=function(v,u,param)
{ rho = param[1];
ga1=param[2]; ga2=param[3]
a1=ga1/sqrt(1-ga1^2)
a2=ga2/sqrt(1-ga2^2)
y1=qsn(u,alpha=a1); y2=qsn(v,alpha=a2)
# parameters on cond cdf to match univariate extended skew-normal
om21=sqrt(1-rho^2)
s=1-ga1^2-ga2^2-rho^2+2*ga1*ga2*rho # should be positive
s=sqrt(s)/om21
tau21=a1*y1
xi21=rho*y1
a21=(ga2-rho*ga1)/s/om21
ccdf=psn(y2,xi=xi21,omega=om21,alpha=a21,tau=tau21)
ccdf
}

#' @note  pcond12 = P(U<u|V=v) for bskewn copula with  y1=qsn(u), y2=qsn(v)
#' @return
#' \code{pcondbskewncop12}: conditional cdf P(U<u|V=v) of bivariate skew-normal copula
#' @rdname bskewncop
#' @import sn CopulaModel
#' @export
pcondbskewncop12=function(u,v,param)
{ rho = param[1];
ga1=param[2]; ga2=param[3]
a1=ga1/sqrt(1-ga1^2)
a2=ga2/sqrt(1-ga2^2)
y1=qsn(u,alpha=a1); y2=qsn(v,alpha=a2)
# parameters on cond cdf to match univariate extended skew-normal
om12=sqrt(1-rho^2)
s=1-ga1^2-ga2^2-rho^2+2*ga1*ga2*rho # should be positive
s=sqrt(s)/om12
tau12=a2*y2
xi12=rho*y2
a12=(ga1-rho*ga2)/s/om12
ccdf=psn(y1,xi=xi12,omega=om12,alpha=a12,tau=tau12)
ccdf
}

#' @note qcond21 = bskewn copula with  y1=qsn(u), y2=qsn(v)
#' @return
#' \code{qcondbskewncop21}: conditional quantile C_{2|1}^{-1}(p|u)of bivariate skew-normal copula
#' @rdname bskewncop
#' @import sn CopulaModel
#' @export
qcondbskewncop21=function(p,u,param)
{ rho = param[1];
ga1=param[2]; ga2=param[3]
a1=ga1/sqrt(1-ga1^2)
a2=ga2/sqrt(1-ga2^2)
y1=qsn(u,alpha=a1); #y2=qsn(v,alpha=a2)
# parameters on cond cdf to match univariate extended skew-normal
om21=sqrt(1-rho^2)
s=1-ga1^2-ga2^2-rho^2+2*ga1*ga2*rho # should be positive
s=sqrt(s)/om21
tau21=a1*y1
xi21=rho*y1
a21=(ga2-rho*ga1)/s/om21
y2=qsn(p,xi=xi21,omega=om21,alpha=a21,tau=tau21)
v=psn(y2,alpha=a2)
v
}

#' @note qcond12 = bskewn copula with  y1=qsn(u), y2=qsn(v)
#' @return
#' \code{qcondbskewncop12}: conditional quantile C_{1|2}^{-1}(p|v)of bivariate skew-normal copula
#' @rdname bskewncop
#' @import sn CopulaModel
#' @export
qcondbskewncop12=function(p,v,param)
{ rho = param[1];
ga1=param[2]; ga2=param[3]
a1=ga1/sqrt(1-ga1^2)
a2=ga2/sqrt(1-ga2^2)
#y1=qsn(u,alpha=a1);
y2=qsn(v,alpha=a2)
# parameters on cond cdf to match univariate extended skew-normal
om12=sqrt(1-rho^2)
s=1-ga1^2-ga2^2-rho^2+2*ga1*ga2*rho # should be positive
s=sqrt(s)/om12
tau12=a2*y2
xi12=rho*y2
a12=(ga1-rho*ga2)/s/om12
y1=qsn(p,xi=xi12,omega=om12,alpha=a12,tau=tau12)
u=psn(y1,alpha=a1)
u
}
