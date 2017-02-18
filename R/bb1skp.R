#' Skew BB1skp Copula Family
#'
#' "p" stands for permutation (u and v are switched from BB1sk).
#' Adds a skew parameter alpha in [0,1]. Distribution function is found by
#' \code{pbb1v(u/u^alpha,v,cpar)*u^alpha}.
#' @author Harry Joe; Vincenzo Coia
#' @rdname bb1skp
#' @export
pbb1skp=function(u,v,cparsk) pbb1sk(v, u, cparsk)


#' @rdname bb1skp
#' @export
pcondbb1skp21=function(v,u,cparsk) pcondbb1sk12(v, u, cparsk)

#' @rdname bb1skp
#' @export
pcondbb1skp <- pcondbb1skp21

#' @rdname bb1skp
#' @export
pcondbb1skp12=function(u,v,cparsk) pcondbb1sk21(u, v, cparsk)

#' @rdname bb1skp
#' @export
dbb1skp=function(u,v,cparsk) dbb1sk(v, u, cparsk)

