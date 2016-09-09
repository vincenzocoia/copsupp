
# modified to get quantile function
# this code is not efficient but should be OK for case of 
#  pair-copulas all permutation symmetric

# Different copula family for each edge
# nsim = #replications or sample size
# parvec = vector of parameters to be optimized in nllk 
# A = dxd vine array with 1:d on diagonal
# ntrunc = truncated level, assume >=1
# pcondmat = matrix of names of conditional cdfs for trees 1,...,ntrunc
#  (assuming only one needed for permutation symmetric pair-copulas)
# qcondmat = matrix of names of conditional quantile functions for 
#        trees 1,...,ntrunc
#   pcondmat and qcondmat are empty for diagonal and lower triangle,
#    and could have ntrunc rows or be empty for rows ntrunc+1 to d-1
# np = dxd where np[ell,j] is #parameters for parameter th[ell,j]
#   for pair-copula in tree ell, variables j and A[ell,j] 
#   np=0 on and below diagonal
# iinv=T to check that this is inverse of rvinepcond()
#    in this case, columns of pmat come from rvinepcond()
# iinv=F, get quantiles C_{d|1:(d-1)}(p|u[1:(d-1)]) based on last column of 
#       pmat[,d] where u[1],..,u[d-1] have been previously converted to
#       u[1], C_{2|1}(u[2]|u[1]), ... C_{d-1|1:(d-2)}(u[d-1]|u[1:(d-2)])
#       via rvinepcond()
## Output: nsim x d matrix with values in (0,1) or
#        quantile C_{d|1...d-1}(p| u1,...,u[d-1])
#rvinesimvec2=function(nsim,A,ntrunc,parvec,np,qcondmat,pcondmat,iprint=F)
temrvineqcond=function(pmat,A,ntrunc,parvec,np,qcondmat,pcondmat,iinv=F)
{ d=ncol(A)
  # get matrix ip1,ip2 of indices
  ii=0
  ip1=matrix(0,d,d); ip2=matrix(0,d,d)
  for(ell in 1:ntrunc)
  { for(j in (ell+1):d)
    { ip1[ell,j]=ii+1; ip2[ell,j]=ii+np[ell,j]
      ii=ii+np[ell,j]
    }
  }
  #if(iprint) { print(ip1); print(ip2) }
  out=varray2M(A)
  M=out$mxarray
  icomp=out$icomp
  #p=matrix(runif(nsim*d),nsim,d)
  p=pmat
  nsim=nrow(p)
  qq=array(0,c(nsim,d,d)); v=array(0,c(nsim,d,d))
  u=matrix(0,nsim,d)
  u[,1]=p[,1]; qq[,1,1]=p[,1]; qq[,2,2]=p[,2];
  qcond=match.fun(qcondmat[1,2])
  u[,2]=qcond(p[,2],p[,1],parvec[ip1[1,2]:ip2[1,2]])
  qq[,1,2]=u[,2]
  if(icomp[1,2]==1) 
  { #pcond=match.fun(pcondnames[1])
    pcond=match.fun(pcondmat[1,2])
    v[,1,2]=pcond(u[,1],u[,2],parvec[ip1[1,2]:ip2[1,2]])
  }
  # the main loop 
  for(j in 3:d)  # variable
  { tt=min(ntrunc,j-1)
    qq[,tt+1,j]=p[,j]
    if(tt>1)
    { for(ell in seq(tt,2)) # tree
      { if(A[ell,j]==M[ell,j]) { s= qq[,ell,A[ell,j]] } 
        else { s=v[,ell-1,M[ell,j]] }
        #qcond=match.fun(qcondnames[ell])
        qcond=match.fun(qcondmat[ell,j])
        qq[,ell,j]=qcond(qq[,ell+1,j], s, parvec[ip1[ell,j]:ip2[ell,j]]); 
      }
    }
    #qcond=match.fun(qcondnames[1])
    qcond=match.fun(qcondmat[1,j])
    qq[,1,j]=qcond(qq[,2,j],u[,A[1,j]], parvec[ip1[1,j]:ip2[1,j]])
    u[,j]=qq[,1,j] 
    # set up for next iteration (not needed for last j=d)
    #pcond=match.fun(pcondnames[1])
    pcond=match.fun(pcondmat[1,j])
    v[,1,j]=pcond(u[,A[1,j]],u[,j],parvec[ip1[1,j]:ip2[1,j]])
    if(tt>1)
    { for(ell in 2:tt)
      { if(A[ell,j]==M[ell,j]) { s=qq[,ell,A[ell,j]] }
        else { s=v[,ell-1,M[ell,j]] }
        if(icomp[ell,j]==1) 
        { #pcond=match.fun(pcondnames[ell])
          pcond=match.fun(pcondmat[ell,j])
          v[,ell,j]=pcond(s, qq[,ell,j], parvec[ip1[ell,j]:ip2[ell,j]]); 
        }
      }
    }
  }
  if(iinv) u=u[,d]
  u
}

#============================================================

library(CopulaModel)
source("hjcondvine.R")
d=5
np=matrix(1,d,d)
ntrunc=4
udat=matrix(c(.1,.2,.3,.4,.5, .6,.7,.8,.9,.5),2,5,byrow=T)
parvec=c(2,2,2,2,1.5,1.5,1.5,1.2,1.2,1.1)
D5=Dvinearray(d)
pcondmat=matrix("pcondgum",d,d)
qcondmat=matrix("qcondgum",d,d)

pmatd=rvinepcond(parvec,udat,D5,ntrunc,pcondmat,np)
print(pmatd)

#rvineqcond=function(pmat,A,ntrunc,parvec,np,qcondmat,pcondmat,iprint=F)

qmatd=temrvineqcond(pmatd,D5,ntrunc,parvec,np,qcondmat,pcondmat,iinv=F)
print(qmatd) # same as udat

q5=temrvineqcond(pmatd,D5,ntrunc,parvec,np,qcondmat,pcondmat,iinv=T)
print(q5)
# [1] 0.5 0.5
pnewd=pmatd; pnewd[,d]=0.9
q5b=temrvineqcond(pnewd,D5,ntrunc,parvec,np,qcondmat,pcondmat,iinv=T)
print(q5b)
# [1] 0.5227173 0.8994989

# direct approach
ntrunc=4
pcond=pcondgum
qcond=qcondgum
parm=matrix(0,d,d) # format of parameter for vine
ii=0;
for(ell in 1:ntrunc)
{ for(j in (ell+1):d)
  { parm[ell,j]=parvec[ii+1]
    ii=ii+1
  }
}

# compared with specific case for dvine
pdchk=udat
pdchk[,2]=pcond(udat[,2],udat[,1],parm[1,2])
tem12=pcond(udat[,1],udat[,2],parm[1,2])
tem32=pcond(udat[,3],udat[,2],parm[1,3])
pdchk[,3]=pcond(tem32,tem12,parm[2,3])
tem1g23=pcond(tem12,tem32,parm[2,3])
tem23=pcond(udat[,2],udat[,3],parm[1,3])
tem43=pcond(udat[,4],udat[,3],parm[1,4])
tem4g23=pcond(tem43,tem23,parm[2,4])
pdchk[,4]=pcond(tem4g23,tem1g23,parm[3,4])
tem1g234=pcond(tem1g23,tem4g23,parm[3,4])
tem34=pcond(udat[,3],udat[,4],parm[1,4])
tem54=pcond(udat[,5],udat[,4],parm[1,5])
tem5g34=pcond(tem54,tem34,parm[2,5])
tem2g34=pcond(tem23,tem43,parm[2,4])
tem5g234=pcond(tem5g34,tem2g34,parm[3,5])
pdchk[,5]=pcond(tem5g234,tem1g234,parm[4,5])
# quantiles at .9
p=rep(.9,2)
s2=qcond(p,tem1g234,parm[4,5])
s3=qcond(s2,tem2g34,parm[3,5])
s4=qcond(s3,tem34,parm[2,5])
s5=qcond(s4,udat[,4],parm[1,5])
print(s5)
# 1] 0.5227173 0.8994989

#============================================================

# truncated vines

# ntrunc=2
pmatd2=rvinepcond(parvec,udat,D5,2,pcondmat,np)
#     [,1]      [,2]      [,3]      [,4]       [,5]
#[1,]  0.1 0.4938008 0.7385935 0.7430672 0.76377507
#[2,]  0.6 0.7328919 0.8830648 0.9508236 0.08757126
qmatd2=temrvineqcond(pmatd2,D5,2,parvec,np,qcondmat,pcondmat,iinv=F)
print(qmatd2) # same as udat
q52=temrvineqcond(pmatd2,D5,2,parvec,np,qcondmat,pcondmat,iinv=T)
print(q52)
# [1] 0.5 0.5
pnewd=pmatd2; pnewd[,d]=0.9
q52b=temrvineqcond(pnewd,D5,2,parvec,np,qcondmat,pcondmat,iinv=T)
print(q52b)
# [1] 0.6190225 0.9281953

# ntrunc=1
pmatd1=rvinepcond(parvec,udat,D5,1,pcondmat,np)
#     [,1]      [,2]      [,3]      [,4]       [,5]
#[1,]  0.1 0.4938008 0.5364857 0.5842195 0.63198274
#[2,]  0.6 0.7328919 0.7951641 0.8831572 0.08282512
qmatd1=temrvineqcond(pmatd1,D5,1,parvec,np,qcondmat,pcondmat,iinv=F)
print(qmatd1) # same as udat
q51=temrvineqcond(pmatd1,D5,1,parvec,np,qcondmat,pcondmat,iinv=T)
print(q51)
# [1] 0.5 0.5
pnewd=pmatd1; pnewd[,d]=0.9
q51b=temrvineqcond(pnewd,D5,1,parvec,np,qcondmat,pcondmat,iinv=T)
print(q51b)
# [1] 0.7332716 0.9529803

