#Use approximated K inverse to estimate loading matrix Z
#
#
#
#
#
########################################################################################################################
# Package: SpatialPCA
# Version: 1.1.0
# Date   : 2021-10-27
# Title : Spatially Aware Dimension Reduction for Spatial Transcriptomics
# Authors: L. Shang and X. Zhou
# Contacts: shanglu@umich.edu
#          University of Michigan, Department of Biostatistics
####################################################################################################


#' Calculate loading matrix.
#'
#' @param object SpatialPCA object.
#' @param maxiter Maximum iteration number. Default is 300.
#' @param initial_tau Initial value of tau. Default is 1. Because we need tau to be positive, we calculate exp(log(tau)) during iterations.
#' @param fast Select "TRUE" if the user wants to use low-rank approximation on the kernel matrix to accelerate the algorithm, otherwise select "FALSE".
#' @param eigenvecnum When fast=TRUE, eigenvecnum is the number of top eigenvectors and eigenvalues to be used in low-rank approximation in the eigen decomposition step for kernel matrix.
#' The default is NULL, if specified, it is recommended to use eigenvecnum=20 when sample size is large (e.g. >5,000). When sample size is small, eigenvecnum is suggested to explain at least 90% variance.
#' @param SpatialPCnum Number of spatial PCs.
#' 
#' @return Returns SpatialPCA object with estimated loading matrix W.
#'
#' @import RSpectra
#'
#' @export
#'
SpatialPCA_EstimateLoading_NNGP = function(object, maxiter=300,initial_tau=1,fast=FALSE,eigenvecnum=NULL,SpatialPCnum=20, sigma_inverse, sigma_inverse_det, sigma_det){
  
  suppressMessages(require(RSpectra))
  set.seed(1234)
  param_ini=log(initial_tau)
  object@SpatialPCnum = SpatialPCnum
  object@fast = fast
  object@params$X = scale(object@location)
  object@params$n = dim(object@params$X)[1]
  object@params$p=dim(object@params$X)[2]
  
  if(is.null(object@covariate)){
    object@params$H = matrix(1, dim(object@params$X)[1],1)
    HH_inv=solve(t(object@params$H)%*%object@params$H,tol = 1e-40)
    HH = object@params$H%*%HH_inv%*%t(object@params$H)
    object@params$M=diag(object@params$n)-HH
    # Y=expr
    object@params$tr_YMY=sum(diag(object@params$expr%*%object@params$M%*%t(object@params$expr)))
    object@params$YM = object@params$expr%*%object@params$M
    object@params$q=1
  }else{
    object@params$q = dim(object@covariate)[2]+1
    object@params$H = matrix(0, object@params$n,object@params$q)
    object@params$H[,1]=1
    object@params$H[,2:object@params$q] = object@covariate
    HH_inv=solve(t(object@params$H)%*%object@params$H,tol = 1e-40)
    HH=object@params$H%*%HH_inv%*%t(object@params$H)
    object@params$M=diag(object@params$n)-HH
    #Y=expr
    object@params$tr_YMY=sum(diag(object@params$expr%*%object@params$M%*%t(object@params$expr)))
    object@params$YM = object@params$expr%*%object@params$M
  }
  
  
  object@params$MYt = object@params$M %*% t(object@params$expr)
  object@params$YMMYt = object@params$YM %*% object@params$MYt
  object@params$Xt = t(object@params$H)
  object@params$SpatialPCnum = SpatialPCnum
  
  
  optim_result =try(optim(param_ini, SpatialPCA_estimate_parameter_NNGP,params=object@params,control = list(maxit = maxiter), lower = -10, upper = 10,method="Brent"),silent=T)
  
  object@tau = exp(optim_result$par)
  k = dim(object@params$expr)[1]
  n = dim(object@params$expr)[2]
  q=object@params$q
 
 
  G_each = object@params$YM %*% solve(object@params$M + (1/object@tau * sigma_inverse)) %*% object@params$MYt
  object@W = eigs_sym(G_each, k=SpatialPCnum, which = "LM")$vectors
  object@sigma2_0 = as.numeric((object@params$tr_YMY+F_funct_sameG(object@W,G_each))/(k*(n-q)))

  return(object)
}


#' @import RSpectra
SpatialPCA_estimate_parameter_NNGP = function(param_ini, params){
  # suppressMessages(require(RSpectra))
  set.seed(1234)
  tau=exp(param_ini[1])
  k = dim(object@params$expr)[1]
  n = dim(object@params$expr)[2]
  q=object@params$q
  PCnum=object@params$SpatialPCnum

  G_each = object@params$YM %*% solve(object@params$M + (1/tau * as(sigma_inverse,"sparseMatrix"))) %*% object@params$MYt
  
  log_det_tauK_I = log(tau*sigma_det + object@params$n*1)
 
  Xt_invmiddle_X = object@params$Xt%*%solve(diag(object@params$n) + tau * object@kernelmat)%*% object@params$H
  
  log_det_Xt_inv_X = determinant(Xt_invmiddle_X, logarithm=TRUE)$modulus[1]
  
  sum_det=0
  
  sum_det=sum_det+(0.5*log_det_tauK_I+0.5*log_det_Xt_inv_X  )*PCnum
  
  W_est_here = eigs_sym(G_each, k=PCnum, which = "LM")$vectors
  -(-sum_det -(k*(n-q))/2*log(object@params$tr_YMY+F_funct_sameG(W_est_here,G_each)))
}


F_funct_sameG = function(X,G){ # G is a matrix
  return_val=0
  for(i in 1: dim(X)[2]){
    return_val=return_val+t(X[,i])%*%G%*%X[,i]
  }
  -return_val
}
