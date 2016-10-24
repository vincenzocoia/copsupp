#' @param cparsk copula parameter vector, third component in [0,1]
#' @param eps tolerance for convergence
#' @param mxiter maximum number of Newton-Raphson iterations
#' @rdname bb1sk
#' @export
qcondbb1sk21=function(p,u,cparsk, eps=1.e-8,mxiter=30,iprint=F)
{
  iter=0; diff=1.;
  v=qcondbb1(p,u,cparsk[1:2])  # starting point
  #while(iter<mxiter & abs(diff)>eps)
  while(iter<mxiter & max(abs(diff))>eps)
  { num=pcondbb1sk21(v,u,cparsk)-p;
    den=dbb1sk(u,v,cparsk);
    diff=num/den;
    v=v-diff;
    #while(v<0. || v>1.) { diff=diff/2.; v=v+diff;}
    while(min(v)<0. | max(v)>1.) { diff=diff/2.; v=v+diff;}
    iter=iter+1;
    if(iprint) cat(iter,v,"\n")
  }
  v
}

#' @rdname bb1sk
#' @export
qcondbb1sk12=function(p,v,cparsk, eps=1.e-8,mxiter=30,iprint=F)
{
  iter=0; diff=1.;
  u=qcondbb1(p,v,cparsk[1:2])  # starting point
  #while(iter<mxiter & abs(diff)>eps)
  while(iter<mxiter & max(abs(diff))>eps)
  { num=pcondbb1sk12(u,v,cparsk)-p;
    den=dbb1sk(u,v,cparsk);
    diff=num/den;
    u=u-diff;
    #while(u<0. || u>1.) { diff=diff/2.; u=u+diff;}
    while(min(u)<0. | max(u)>1.) { diff=diff/2.; u=u+diff;}
    iter=iter+1;
    if(iprint) cat(iter,u,"\n")
  }
  u
}

#' @rdname bb1sk
#' @export
qcondbb1sk <- qcondbb1sk21
