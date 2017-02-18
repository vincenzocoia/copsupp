#' @param cparsk copula parameter vector, third component in [0,1]
#' @param eps tolerance for convergence
#' @param mxiter maximum number of Newton-Raphson iterations
#' @rdname bb1skp
#' @export
qcondbb1skp21=function(p,u,cparsk, eps=1.e-8,mxiter=30,iprint=F)
    qcondbb1sk12(p,u,cparsk, eps=eps,mxiter=mxiter,iprint=iprint)

#' @rdname bb1skp
#' @export
qcondbb1skp12=function(p,v,cparsk, eps=1.e-8,mxiter=30,iprint=F)
    qcondbb1sk21(p,u,cparsk, eps=eps,mxiter=mxiter,iprint=iprint)

#' @rdname bb1skp
#' @export
qcondbb1skp <- qcondbb1skp21
