#' # Assume vine array with 1:d on diagonal.
#' # Vine functions for mapping of (u[1],...,u[d])  to
#' # p[1]=u[1], p[2]=C_{2|1}(u[2]|u[1]) , ...
#' #         p[d]=C_{d|1:(d-1)}(u[d]}u[1:(d-1)])
#' # and inverse.
#' # Byproducts :
#' #  (a) conditional cdf F_{d|1:(d-1)}(x_d| x_1,...,x_{d-1}) after converting
#' #   x_j to u_j=F_j(x_j) for j=1,...,d
#' #  (b) conditional quantile F^{-1}_{d|1:(d-1)}(p|x_1,...,x_{d-1})
#' #  after converting x_j to u_jF_j(x_j) for j=1,...,d-1 and
#' #  applying F_d^{-1}(u_d) to the output u_d
#'
#' # Note that if pair-copula families are permutation asymmetric,
#' # the code below needs to replace pcond (pcondnames)
#' # with pcond21 (pcond21names) and pcond12 (pcond12names),
#' # with small changes to the code, replacing pcond by either pcond21 or pcond12
#'
#' # so below is not valid for pcondgumv pcondgumu etc where
#' # one of the U(0,1) variables has been reflected to ger negative depc.
#'
#' #============================================================
#'
#' # Versions of pcond vectors for R-vine
#' #   with different pair-copula family for each edge of the vine
#' # for each row u[1:d] of hypercube data,
#' # return
#' # u1, C_{2|1}(u2|u1), C_{3|12}(u3|u1,u2), ...  C_{d|1..d-1}(u_d|u[1:(d-1)])
#'
#' #' Copula Conditional Distributions
#' #'
#' #' These functions compute the conditional distributions of the
#' #' last variable, given the other variables in the vine.
#' #'
#' #' parvec = vector of parameters of pair-copulas
#' #' udat = nxd matrix with uniform scores
#' #' A = dxd vine array with 1:d on diagonal
#' #' ntrunc = truncated level, assume >=1
#' #' pcondmat = matrix of names of conditional cdfs for trees 1,...,ntrunc
#' #'  (assuming only one needed for permutation symmetric pair-copulas)
#' #'   pcondmat is empty for diagonal and lower triangle,
#' #'    and could have ntrunc rows or be empty for rows ntrunc+1 to d-1
#' #' np = dxd where np[ell,j] is size for parameter th[ell,j]
#' #'   for pair-copula in tree ell, variables j and A[ell,j]
#' #'   np=0 on and below diagonal
#' #' Output:
#' #' return C_{2|1}(u2|u1), C_{3|12}(u3|u1,u2), ...  C_{d|1..d-1}(ud|u1,...u[d-1])
#' #' [This function is modification of rvinellkv.trunc2 in CopulaModel]
#' #' @author Harry Joe
#' #' @import CopulaModel
#' #' @export
#' #rvinellkv.trunc2=function(parvec,udat,A,ntrunc,logdcopmat,pcondmat,np)
#' rvinepcond=function(parvec,udat,A,ntrunc,pcondmat,np)
#' { d=ncol(A)  # or ncol(udat)
#'   npmax=max(np); th=array(0,c(npmax,d,d)) # format of parameter for vine
#'   ii=0;
#'   for(ell in 1:ntrunc)
#'   { for(j in (ell+1):d)
#'     { ipp=1:np[ell,j]
#'       th[ipp,ell,j]=parvec[ii+ipp]
#'       ii=ii+np[ell,j]
#'     }
#'   }
#'   out=varray2M(A)
#'   M=out$mxarray
#'   icomp=out$icomp
#'   n=nrow(udat)
#'   #llkv=rep(0,n)
#'   vinepcond=matrix(0,n,d)
#'   vinepcond[,1]=udat[,1]
#'   v=matrix(0,n,d); vp=matrix(0,n,d); s=matrix(0,n,d);
#'   nllk=0
#'   # tree 1
#'   #for(j in 2:d)
#'   #{ ipp=1:np[1,j]
#'     #logdcop=match.fun(logdcopmat[1,j])
#'     #llkv=llkv+logdcop(udat[,A[1,j]],udat[,j],th[ipp,1,j])
#'   #}
#'   # tree 2
#'   if(ntrunc>=2)
#'   { for(j in 2:d)
#'     { ipp=1:np[1,j]
#'       pcond=match.fun(pcondmat[1,j])
#'       if(icomp[1,j]==1) vp[,j]=pcond(udat[,A[1,j]],udat[,j],th[ipp,1,j])
#'       v[,j]=pcond(udat[,j],udat[,A[1,j]],th[ipp,1,j])
#'     }
#'     vinepcond[,2]=v[,2]
#'     for(j in 3:d)
#'     { if(A[2,j]<M[2,j]) s[,j]=vp[,M[2,j]] else s[,j]=v[,A[2,j]] }
#'     #for(j in 3:d)
#'     #{ ipp=1:np[2,j]
#'       #logdcop=match.fun(logdcopmat[2,j])
#'       #llkv=llkv+logdcop(s[,j],v[,j],th[ipp,2,j])
#'     #}
#'     w=v; wp=vp
#'   }
#'   # remaining trees
#'   if(ntrunc>=3)
#'   { for(ell in 3:ntrunc)
#'     { for(j in ell:d)
#'       { ipp=1:np[ell-1,j]
#'         pcond=match.fun(pcondmat[ell-1,j])
#'         if(icomp[ell-1,j]==1) vp[,j]=pcond(s[,j],w[,j],th[ipp,ell-1,j])
#'         v[,j]=pcond(w[,j],s[,j],th[ipp,ell-1,j])
#'       }
#'       vinepcond[,ell]=v[,ell]
#'       for(j in (ell+1):d)
#'       { if(A[ell,j]<M[ell,j]) s[,j]=vp[,M[ell,j]] else s[,j]=v[,A[ell,j]] }
#'       #for(j in (ell+1):d)
#'       #{ ipp=1:np[ell,j]
#'         #logdcop=match.fun(logdcopmat[ell,j])
#'         #llkv=llkv+logdcop(s[,j],v[,j],th[ipp,ell,j])
#'       #}
#'       w=v; # wp=vp  # wp is not used
#'     }
#'   }
#'   if(ntrunc>1 & ntrunc<d)
#'   { ell=ntrunc+1
#'     for(j in ell:d)
#'     { ipp=1:np[ell-1,j]
#'       pcond=match.fun(pcondmat[ell-1,j])
#'       if(icomp[ell-1,j]==1) vp[,j]=pcond(s[,j],w[,j],th[ipp,ell-1,j])
#'       v[,j]=pcond(w[,j],s[,j],th[ipp,ell-1,j])
#'       vinepcond[,j]=v[,j]
#'     }
#'   }
#'   if(ntrunc==1)
#'   { ell=2
#'     for(j in ell:d)
#'     { ipp=1:np[ell-1,j]
#'       pcond=match.fun(pcondmat[ell-1,j])
#'       v[,j]=pcond(udat[,j],udat[,A[1,j]],th[ipp,ell-1,j])
#'       vinepcond[,j]=v[,j]
#'     }
#'   }
#'   #llkv
#'   vinepcond
#' }
#'
#' #============================================================
#'
#' # modification of previous R-vine simualtion code to get quantile function
#' # this code is not efficient but should be OK for case of
#' #  pair-copulas all permutation symmetric
#'
#' #' n = #replications or sample size
#' #' pmat = nxd matrix of levels of u1, C_{2|1}(u_2|u_1),...,
#' #'                          C_{d|1..d-1}(u_d|u[1:(d-1)])
#' #' parvec = vector of parameters to be optimized in nllk
#' #' A = dxd vine array with 1:d on diagonal
#' #' ntrunc = truncated level, assume >=1
#' #' pcondmat = matrix of names of conditional cdfs for trees 1,...,ntrunc
#' #'  (assuming only one needed for permutation symmetric pair-copulas)
#' #' qcondmat = matrix of names of conditional quantile functions for
#' #'        trees 1,...,ntrunc
#' #'   qcondmat and pcondmat are empty for diagonal and lower triangle,
#' #'    and could have ntrunc rows or be empty for rows ntrunc+1 to d-1
#' #' np = dxd where np[ell,j] is #parameters for parameter th[ell,j]
#' #'   for pair-copula in tree ell, variables j and A[ell,j]
#' #'   np=0 on and below diagonal
#' #' iinv=T to check that this is inverse of rvinepcond()
#' #'    in this case, columns of pmat come from rvinepcond()
#' #' iinv=F, get quantiles C_{d|1:(d-1)}(p|u[1:(d-1)]) based on last column of
#' #'       pmat[,d] where u[1],..,u[d-1] have been previously converted to
#' #'       u[1], C_{2|1}(u[2]|u[1]), ... C_{d-1|1:(d-2)}(u[d-1]|u[1:(d-2)])
#' #'       via rvinepcond()
#' #' @return n x d matrix with values in (0,1) or
#' #'       quantile C_{d|1...d-1}(p| u[1:(d-1)])
#' #' [This function is modification of rvinesimvec2 in CopulaModel]
#' #' @author Harry Joe
#' #' @import CopulaModel
#' #' @export
#' #rvinesimvec2=function(nsim,A,ntrunc,parvec,np,qcondmat,pcondmat,iprint=F)
#' rvineqcond=function(pmat,A,ntrunc,parvec,np,qcondmat,pcondmat,iinv=F)
#' { d=ncol(A)
#'   # get matrix ip1,ip2 of indices
#'   ii=0
#'   ip1=matrix(0,d,d); ip2=matrix(0,d,d)
#'   for(ell in 1:ntrunc)
#'   { for(j in (ell+1):d)
#'     { ip1[ell,j]=ii+1; ip2[ell,j]=ii+np[ell,j]
#'       ii=ii+np[ell,j]
#'     }
#'   }
#'   #if(iprint) { print(ip1); print(ip2) }
#'   out=varray2M(A)
#'   M=out$mxarray
#'   icomp=out$icomp
#'   #p=matrix(runif(nsim*d),nsim,d)
#'   p=pmat
#'   n=nrow(p)
#'   qq=array(0,c(n,d,d)); v=array(0,c(n,d,d))
#'   u=matrix(0,n,d)
#'   u[,1]=p[,1]; qq[,1,1]=p[,1]; qq[,2,2]=p[,2];
#'   qcond=match.fun(qcondmat[1,2])
#'   u[,2]=qcond(p[,2],p[,1],parvec[ip1[1,2]:ip2[1,2]])
#'   qq[,1,2]=u[,2]
#'   if(icomp[1,2]==1)
#'   { #pcond=match.fun(pcondnames[1])
#'     pcond=match.fun(pcondmat[1,2])
#'     v[,1,2]=pcond(u[,1],u[,2],parvec[ip1[1,2]:ip2[1,2]])
#'   }
#'   # the main loop
#'   for(j in 3:d)  # variable
#'   { tt=min(ntrunc,j-1)
#'     qq[,tt+1,j]=p[,j]
#'     if(tt>1)
#'     { for(ell in seq(tt,2)) # tree
#'       { if(A[ell,j]==M[ell,j]) { s= qq[,ell,A[ell,j]] }
#'         else { s=v[,ell-1,M[ell,j]] }
#'         #qcond=match.fun(qcondnames[ell])
#'         qcond=match.fun(qcondmat[ell,j])
#'         qq[,ell,j]=qcond(qq[,ell+1,j], s, parvec[ip1[ell,j]:ip2[ell,j]]);
#'       }
#'     }
#'     #qcond=match.fun(qcondnames[1])
#'     qcond=match.fun(qcondmat[1,j])
#'     qq[,1,j]=qcond(qq[,2,j],u[,A[1,j]], parvec[ip1[1,j]:ip2[1,j]])
#'     u[,j]=qq[,1,j]
#'     # set up for next iteration (not needed for last j=d)
#'     #pcond=match.fun(pcondnames[1])
#'     pcond=match.fun(pcondmat[1,j])
#'     v[,1,j]=pcond(u[,A[1,j]],u[,j],parvec[ip1[1,j]:ip2[1,j]])
#'     if(tt>1)
#'     { for(ell in 2:tt)
#'       { if(A[ell,j]==M[ell,j]) { s=qq[,ell,A[ell,j]] }
#'         else { s=v[,ell-1,M[ell,j]] }
#'         if(icomp[ell,j]==1)
#'         { #pcond=match.fun(pcondnames[ell])
#'           pcond=match.fun(pcondmat[ell,j])
#'           v[,ell,j]=pcond(s, qq[,ell,j], parvec[ip1[ell,j]:ip2[ell,j]]);
#'         }
#'       }
#'     }
#'   }
#'   if(iinv) u=u[,d]  # last column as quantiles
#'   u  # nxd matrix to check that this funciton is inverse of rvinepcond
#' }
#'
