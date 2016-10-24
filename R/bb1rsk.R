#' Skew BB1r Copula Family
#'
#' Adds a skew parameter alpha in [0,1]. Distribution function is found by
#' \code{pbb1r(u/u^alpha,v,cpar)*u^alpha}.
#' @author Harry Joe
#' @rdname bb1rsk
#' @export
pbb1rsk=function(u,v,cparsk)
{ alp=cparsk[3]
  cpar=cparsk[1:2]
  ua=u^alp
  cdf=pbb1r(u/ua,v,cpar)*ua
  cdf
}

#' @rdname bb1rsk
#' @export
pcondbb1rsk21=function(v,u,cparsk)
{ alp=cparsk[3]
  cpar=cparsk[1:2]
  ua1=u^(1-alp)
  condcdf=(1-alp)*pcondbb1r(v,ua1,cpar) + alp*pbb1r(ua1,v,cpar)/ua1
  condcdf
}

#' @rdname bb1rsk
#' @export
pcondbb1rsk <- pcondbb1rsk21

#' @rdname bb1rsk
#' @export
pcondbb1rsk12=function(u,v,cparsk)
{ alp=cparsk[3]
  cpar=cparsk[1:2]
  ua=u^alp
  condcdf=pcondbb1r(u/ua,v,cpar)*ua
  condcdf
}

#' @rdname bb1rsk
#' @export
dbb1rsk=function(u,v,cparsk)
{ alp=cparsk[3]
  cpar=cparsk[1:2]
  ua1=u^(1-alp)
  pdf=(1-alp)*dbb1r(ua1,v,cpar) + alp*pcondbb1r(ua1,v,cpar)/ua1
  pdf
}


