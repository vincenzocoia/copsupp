# bivariate skew-normal copula with parameters rho, ga1, ga2
# such that determinant with these 3 correlation is positive definite

# Rmat= [1 ga1 ga2; ga1 1 rho; ga2 rho 1 ] 
# is correlation matrix of (Z0,Z1,Z2)
# the bivariate skew normal distribution is (Y1,Y2) = [Z1,Z2|Z0>0]
# the conditional distribution is univariate extended skew-normal
# with 4 parameters: xi=location, omega=scale, alpha=ga/sqrt(1-ga^2),
# tau=(negative threshold) 

# these functions make use of functions in library(sn) and
# library(CopulaModel)

# param=c(rho,ga1,ga2), constraint is that
# correlation matrix with these parameters should be positive definite
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

# skew-normal copula from Azzalini-DallaValle's bivariate skew-normal
# u,v = values or vectors in (0,1)
# param = (rho,ga1,ga2), with -1<rho<1, -1<ga1<1, -1<ga2<1 
#  and a constraint 1-ga1^2-ga2^2-rho^2+2*ga1*ga2*rho>0

# Output: bivariate skew-normal copula density
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

# Output: bivariate skew-normal copula cdf
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

# pcond21 = P(V<v|U=u) for bskewn copula with  y1=qsn(u), y2=qsn(v)
# Output: conditional cdf P(V<v|U=u) of bivariate skew-normal copula 
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

# pcond12 = P(U<u|V=v) for bskewn copula with  y1=qsn(u), y2=qsn(v)
# Output: conditional cdf P(U<u|V=v) of bivariate skew-normal copula 
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

# qcond21 = bskewn copula with  y1=qsn(u), y2=qsn(v)
# Output: conditional quantile C_{2|1}^{-1}(p|u)of bivariate skew-normal copula 
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

# qcond12 = bskewn copula with  y1=qsn(u), y2=qsn(v)
# Output: conditional quantile C_{1|2}^{-1}(p|v)of bivariate skew-normal copula 
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

#============================================================

# conditional u given v
chkcopderiv2=function(u,vvec,cpar,bcdf,pcond,bpdf,str=" ",eps=1.e-4)
{ nn=length(vvec)
  cat("\n",str, " with parameter ", cpar,"\n")
  cat("u   v    cdf      numccdf     ccdf\n")
  cat("u   v    ccdf     numpdf      pdf\n")
  for(i in 1:nn)
  { v=vvec[i]
    cdf=bcdf(u,v,cpar)
    cdf2=bcdf(u,v+eps,cpar)
    ccdf=pcond(u,v,cpar)
    #cat(u,v,cdf,(cdf2-cdf)/eps,ccdf,"\n")
    tem=(cdf2-cdf)/eps
    cat(u,v,cdf,tem,ccdf,"\n")
    if(abs(tem-ccdf)>eps*30) cat("*** ccdf roundoff???\n")
    ccdf2=pcond(u+eps,v,cpar)
    pdf=bpdf(u,v,cpar)
    #cat(u,v,ccdf,(ccdf2-ccdf)/eps,pdf,"\n")
    tem=(ccdf2-ccdf)/eps
    cat(u,v,ccdf,tem,pdf,"\n")
    if(abs(tem-pdf)>eps*30) cat("*** dcop roundoff???\n")
  }
  invisible(0)
}

