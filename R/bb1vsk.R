#' Skew BB1v Copula Family
#'
#' Adds a skew parameter alpha in [0,1]. Distribution function is found by
#' \code{pbb1v(u/u^alpha,v,cpar)*u^alpha}.
#' @author Harry Joe; Vincenzo Coia
#' @rdname bb1vsk
#' @export
pbb1vsk=function(u,v,cparsk)
{ alp=cparsk[3]
  cpar=cparsk[1:2]
  ua=u^alp
  cdf=pbb1v(u/ua,v,cpar)*ua
  cdf
}

#' @rdname bb1vsk
#' @export
pcondbb1vsk21=function(v,u,cparsk)
{ alp=cparsk[3]
  cpar=cparsk[1:2]
  ua1=u^(1-alp)
  condcdf=(1-alp)*pcondbb1v(v,ua1,cpar) + alp*pbb1v(ua1,v,cpar)/ua1
  ## CopulaModel::pcondbb1 seems to produce values outside of (0,1). For example,
  ##  try CopulaModel::pcondbb1(0.99, 0.99, c(1.8, 0.7))  (=1.28...).
  pmax(condcdf, 0.000001)
}

#' @rdname bb1vsk
#' @export
pcondbb1vsk <- pcondbb1vsk21

#' @rdname bb1vsk
#' @export
pcondbb1vsk12=function(u,v,cparsk)
{ alp=cparsk[3]
  cpar=cparsk[1:2]
  ua=u^alp
  condcdf=pcondbb1v(u/ua,v,cpar)*ua
  condcdf
}

#' @rdname bb1vsk
#' @export
dbb1vsk=function(u,v,cparsk)
{ alp=cparsk[3]
  cpar=cparsk[1:2]
  ua1=u^(1-alp)
  pdf=(1-alp)*dbb1v(ua1,v,cpar) + alp*pcondbb1v(ua1,v,cpar)/ua1
  pdf
}


