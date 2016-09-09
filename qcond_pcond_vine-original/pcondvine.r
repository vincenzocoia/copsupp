
# Versions of pcond vectors for R-vine
#   with different pair-copula family for each edge of the vine
# for each row u[1:d] of hypercube data,
# return C_{2|1}(u2|u1), C_{3|12}(u3|u1,u2), ...  C_{d|1..d-1}(ud|u1,...u[d-1])

# parvec = vector of parameters of pair-copulas
# udat = nxd matrix with uniform scores
# A = dxd vine array with 1:d on diagonal
# ntrunc = truncated level, assume >=1
# pcondmat = matrix of names of conditional cdfs for trees 1,...,ntrunc
#  (assuming only one needed for permutation symmetric pair-copulas)
#   pcondmat is empty for diagonal and lower triangle,
#    and could have ntrunc rows or be empty for rows ntrunc+1 to d-1
# np = dxd where np[ell,j] is size for parameter th[ell,j]
#   for pair-copula in tree ell, variables j and A[ell,j] 
#   np=0 on and below diagonal
# Output: 
# return C_{2|1}(u2|u1), C_{3|12}(u3|u1,u2), ...  C_{d|1..d-1}(u[d]|u[1:(d-1)])
#rvinellkv.trunc2=function(parvec,udat,A,ntrunc,logdcopmat,pcondmat,np)
rvinepcond.trunc=function(parvec,udat,A,ntrunc,pcondmat,np)
{ d=ncol(A)  # or ncol(udat)
  npmax=max(np); th=array(0,c(npmax,d,d)) # format of parameter for vine
  ii=0;
  for(ell in 1:ntrunc)
  { for(j in (ell+1):d)
    { ipp=1:np[ell,j]
      th[ipp,ell,j]=parvec[ii+ipp]
      ii=ii+np[ell,j]
    }
  }
  out=varray2M(A)
  M=out$mxarray
  icomp=out$icomp
  n=nrow(udat)
  #llkv=rep(0,n)
  vinepcond=matrix(0,n,d)
  vinepcond[,1]=udat[,1]
  v=matrix(0,n,d); vp=matrix(0,n,d); s=matrix(0,n,d);
  nllk=0
  # tree 1
  #for(j in 2:d) 
  #{ ipp=1:np[1,j]
    #logdcop=match.fun(logdcopmat[1,j])
    #llkv=llkv+logdcop(udat[,A[1,j]],udat[,j],th[ipp,1,j])
  #}
  # tree 2
  if(ntrunc>=2)
  { for(j in 2:d) 
    { ipp=1:np[1,j]
      pcond=match.fun(pcondmat[1,j])
      if(icomp[1,j]==1) vp[,j]=pcond(udat[,A[1,j]],udat[,j],th[ipp,1,j]) 
      v[,j]=pcond(udat[,j],udat[,A[1,j]],th[ipp,1,j])
    }
    vinepcond[,2]=v[,2]
    for(j in 3:d) 
    { if(A[2,j]<M[2,j]) s[,j]=vp[,M[2,j]] else s[,j]=v[,A[2,j]] } 
    #for(j in 3:d) 
    #{ ipp=1:np[2,j]
      #logdcop=match.fun(logdcopmat[2,j])
      #llkv=llkv+logdcop(s[,j],v[,j],th[ipp,2,j])
    #}
    w=v; wp=vp
  }
  # remaining trees
  if(ntrunc>=3)
  { for(ell in 3:ntrunc)
    { for(j in ell:d) 
      { ipp=1:np[ell-1,j] 
        pcond=match.fun(pcondmat[ell-1,j]) 
        if(icomp[ell-1,j]==1) vp[,j]=pcond(s[,j],w[,j],th[ipp,ell-1,j]) 
        v[,j]=pcond(w[,j],s[,j],th[ipp,ell-1,j])
      }
      vinepcond[,ell]=v[,ell]
      for(j in (ell+1):d) 
      { if(A[ell,j]<M[ell,j]) s[,j]=vp[,M[ell,j]] else s[,j]=v[,A[ell,j]] } 
      #for(j in (ell+1):d) 
      #{ ipp=1:np[ell,j] 
        #logdcop=match.fun(logdcopmat[ell,j]) 
        #llkv=llkv+logdcop(s[,j],v[,j],th[ipp,ell,j])
      #}
      w=v; # wp=vp  # wp is not used 
    }
  }
  if(ntrunc>1 & ntrunc<d)
  { ell=ntrunc+1
    for(j in ell:d) 
    { ipp=1:np[ell-1,j] 
      pcond=match.fun(pcondmat[ell-1,j]) 
      if(icomp[ell-1,j]==1) vp[,j]=pcond(s[,j],w[,j],th[ipp,ell-1,j]) 
      v[,j]=pcond(w[,j],s[,j],th[ipp,ell-1,j])
      vinepcond[,j]=v[,j] 
    }
  }
  if(ntrunc==1)
  { ell=2
    for(j in ell:d) 
    { ipp=1:np[ell-1,j] 
      pcond=match.fun(pcondmat[ell-1,j]) 
      v[,j]=pcond(udat[,j],udat[,A[1,j]],th[ipp,ell-1,j])
      vinepcond[,j]=v[,j] 
    }
  }
  #llkv
  vinepcond
}

#============================================================

library(CopulaModel)
d=5
np=matrix(1,d,d)
ntrunc=4
udat=matrix(c(.1,.2,.3,.4,.5, .6,.7,.8,.9,.5),2,5,byrow=T)
parvec=c(2,2,2,2,1.5,1.5,1.5,1.2,1.2,1.1)
D5=Dvinearray(d)
pcondmat=matrix("pcondgum",d,d)

outd=rvinepcond.trunc(parvec,udat,D5,ntrunc,pcondmat,np)
print(outd)
#     [,1]      [,2]      [,3]      [,4]      [,5]
#[1,]  0.1 0.4938008 0.7385935 0.8312337 0.8782010
#[2,]  0.6 0.7328919 0.8830648 0.9754268 0.1403602

outd2=rvinepcond.trunc(parvec,udat,D5,2,pcondmat,np)
print(outd2)
#     [,1]      [,2]      [,3]      [,4]       [,5]
#[1,]  0.1 0.4938008 0.7385935 0.7430672 0.76377507
#[2,]  0.6 0.7328919 0.8830648 0.9508236 0.08757126

outd1=rvinepcond.trunc(parvec,udat,D5,1,pcondmat,np)
print(outd1)
#     [,1]      [,2]      [,3]      [,4]       [,5]
#[1,]  0.1 0.4938008 0.5364857 0.5842195 0.63198274
#[2,]  0.6 0.7328919 0.7951641 0.8831572 0.08282512


# C-vine
C5=Cvinearray(d)
outc=rvinepcond.trunc(parvec,udat,C5,ntrunc,pcondmat,np)
print(outc)
#     [,1]      [,2]      [,3]      [,4]      [,5]
#[1,]  0.1 0.4938008 0.7228406 0.8433683 0.9167192
#[2,]  0.6 0.7328919 0.8700505 0.9776872 0.1262401

# different truncation levels
outc3=rvinepcond.trunc(parvec,udat,C5,3,pcondmat,np)
print(outc3)
#     [,1]      [,2]      [,3]      [,4]      [,5]
#[1,]  0.1 0.4938008 0.7228406 0.8433683 0.9319560
#[2,]  0.6 0.7328919 0.8700505 0.9776872 0.1922732
outc2=rvinepcond.trunc(parvec,udat,C5,2,pcondmat,np)
print(outc2)
#     [,1]      [,2]      [,3]      [,4]      [,5]
#[1,]  0.1 0.4938008 0.7228406 0.8530008 0.9294355
#[2,]  0.6 0.7328919 0.8700505 0.9828434 0.2832124

outc1=rvinepcond.trunc(parvec,udat,C5,1,pcondmat,np)
print(outc1)
#     [,1]      [,2]      [,3]      [,4]      [,5]
#[1,]  0.1 0.4938008 0.6592740 0.7794878 0.8646421
#[2,]  0.6 0.7328919 0.8746494 0.9689105 0.4179746


#============================================================

parm=matrix(0,d,d) # format of parameter for vine
ii=0;
for(ell in 1:ntrunc)
{ for(j in (ell+1):d)
  { parm[ell,j]=parvec[ii+1]
    ii=ii+1
  }
}

pcond=pcondgum

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
print(pdchk)
#     [,1]      [,2]      [,3]      [,4]      [,5]
#[1,]  0.1 0.4938008 0.7385935 0.8312337 0.8782010
#[2,]  0.6 0.7328919 0.8830648 0.9754268 0.1403602

# truncated vine:
# 1-truncated
pdchk1=udat
pdchk1[,2]=pcond(udat[,2],udat[,1],parm[1,2])
pdchk1[,3]=pcond(udat[,3],udat[,2],parm[1,3])
pdchk1[,4]=pcond(udat[,4],udat[,3],parm[1,4])
pdchk1[,5]=pcond(udat[,5],udat[,4],parm[1,5])
print(pdchk1)
#     [,1]      [,2]      [,3]      [,4]       [,5]
#[1,]  0.1 0.4938008 0.5364857 0.5842195 0.63198274
#[2,]  0.6 0.7328919 0.7951641 0.8831572 0.08282512

# 2-truncated
pdchk2=udat
pdchk2[,2]=pcond(udat[,2],udat[,1],parm[1,2])
tem12=pcond(udat[,1],udat[,2],parm[1,2])
tem32=pcond(udat[,3],udat[,2],parm[1,3])
pdchk2[,3]=pcond(tem32,tem12,parm[2,3])
tem1g23=pcond(tem12,tem32,parm[2,3])
tem23=pcond(udat[,2],udat[,3],parm[1,3])
tem43=pcond(udat[,4],udat[,3],parm[1,4])
tem4g23=pcond(tem43,tem23,parm[2,4])
pdchk2[,4]=tem4g23
tem1g234=tem1g23
tem34=pcond(udat[,3],udat[,4],parm[1,4])
tem54=pcond(udat[,5],udat[,4],parm[1,5])
tem5g34=pcond(tem54,tem34,parm[2,5])
tem2g34=pcond(tem23,tem43,parm[2,4])
tem5g234=tem5g34
pdchk2[,5]=tem5g234
print(pdchk2)
#     [,1]      [,2]      [,3]      [,4]       [,5]
#[1,]  0.1 0.4938008 0.7385935 0.7430672 0.76377507
#[2,]  0.6 0.7328919 0.8830648 0.9508236 0.08757126

#
# compared with specific case for cvine

pcchk=udat
tem21=pcond(udat[,2],udat[,1],parm[1,3])
pcchk[,2]=tem21
tem31=pcond(udat[,3],udat[,1],parm[1,3])
tem3g12=pcond(tem31,tem21,parm[2,3])
pcchk[,3]=tem3g12
tem41=pcond(udat[,4],udat[,1],parm[1,4])
tem4g12=pcond(tem41,tem21,parm[2,4])
tem4g123=pcond(tem4g12,tem3g12,parm[3,4])
pcchk[,4]=tem4g123
tem51=pcond(udat[,5],udat[,1],parm[1,5])
tem5g12=pcond(tem51,tem21,parm[2,5])
tem5g123=pcond(tem5g12,tem3g12,parm[3,5])
pcchk[,5]=pcond(tem5g123,tem4g123,parm[4,5])

print(pcchk)
#     [,1]      [,2]      [,3]      [,4]      [,5]
#[1,]  0.1 0.4938008 0.7228406 0.8433683 0.9167192
#[2,]  0.6 0.7328919 0.8700505 0.9776872 0.1262401

# truncated cases

# 3-truncated
pcchk3=udat
tem21=pcond(udat[,2],udat[,1],parm[1,3])
pcchk3[,2]=tem21
tem31=pcond(udat[,3],udat[,1],parm[1,3])
tem3g12=pcond(tem31,tem21,parm[2,3])
pcchk3[,3]=tem3g12
tem41=pcond(udat[,4],udat[,1],parm[1,4])
tem4g12=pcond(tem41,tem21,parm[2,4])
tem4g123=pcond(tem4g12,tem3g12,parm[3,4])
pcchk3[,4]=tem4g123
tem51=pcond(udat[,5],udat[,1],parm[1,5])
tem5g12=pcond(tem51,tem21,parm[2,5])
tem5g123=pcond(tem5g12,tem3g12,parm[3,5])
pcchk3[,5]=tem5g123
print(pcchk3)
#     [,1]      [,2]      [,3]      [,4]      [,5]
#[1,]  0.1 0.4938008 0.7228406 0.8433683 0.9319560
#[2,]  0.6 0.7328919 0.8700505 0.9776872 0.1922732

# 2-truncated
pcchk2=udat
tem21=pcond(udat[,2],udat[,1],parm[1,3])
pcchk2[,2]=tem21
tem31=pcond(udat[,3],udat[,1],parm[1,3])
tem3g12=pcond(tem31,tem21,parm[2,3])
pcchk2[,3]=tem3g12
tem41=pcond(udat[,4],udat[,1],parm[1,4])
tem4g12=pcond(tem41,tem21,parm[2,4])
tem4g123=tem4g12
pcchk2[,4]=tem4g123
tem51=pcond(udat[,5],udat[,1],parm[1,5])
tem5g12=pcond(tem51,tem21,parm[2,5])
tem5g123=tem5g12
pcchk2[,5]=tem5g123
print(pcchk2)
#     [,1]      [,2]      [,3]      [,4]      [,5]
#[1,]  0.1 0.4938008 0.7228406 0.8530008 0.9294355
#[2,]  0.6 0.7328919 0.8700505 0.9828434 0.2832124

# 1-truncated
pcchk1=udat
pcchk1[,2]=pcond(udat[,2],udat[,1],parm[1,2])
pcchk1[,3]=pcond(udat[,3],udat[,1],parm[1,3])
pcchk1[,4]=pcond(udat[,4],udat[,1],parm[1,4])
pcchk1[,5]=pcond(udat[,5],udat[,1],parm[1,5])
print(pcchk1)
#     [,1]      [,2]      [,3]      [,4]      [,5]
#[1,]  0.1 0.4938008 0.6592740 0.7794878 0.8646421
#[2,]  0.6 0.7328919 0.8746494 0.9689105 0.4179746

