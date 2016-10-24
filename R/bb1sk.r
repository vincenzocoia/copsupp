#' Skew BB1 Copula Family
#'
#' Adds a skew parameter alpha in [0,1]. Distribution function is found by
#' \code{pbb1(u/u^alpha,v,cpar)*u^alpha}.
#' @author Vincenzo Coia, Harry Joe
#' @rdname bb1sk
#' @export
pbb1sk=function(u,v,cparsk)
{ alp=cparsk[3]
  cpar=cparsk[1:2]
  ua=u^alp
  cdf=pbb1(u/ua,v,cpar)*ua
  cdf
}

#' @rdname bb1sk
#' @export
pcondbb1sk21=function(v,u,cparsk)
{ alp=cparsk[3]
  cpar=cparsk[1:2]
  ua1=u^(1-alp)
  condcdf=(1-alp)*pcondbb1(v,ua1,cpar) + alp*pbb1(ua1,v,cpar)/ua1
  condcdf
}

#' @rdname bb1sk
#' @export
pcondbb1sk <- pcondbb1sk21

#' @rdname bb1sk
#' @export
pcondbb1sk12=function(u,v,cparsk)
{ alp=cparsk[3]
  cpar=cparsk[1:2]
  ua=u^alp
  condcdf=pcondbb1(u/ua,v,cpar)*ua
  condcdf
}

#' @rdname bb1sk
#' @export
dbb1sk=function(u,v,cparsk)
{ alp=cparsk[3]
  cpar=cparsk[1:2]
  ua1=u^(1-alp)
  pdf=(1-alp)*dbb1(ua1,v,cpar) + alp*pcondbb1(ua1,v,cpar)/ua1
  pdf
}


