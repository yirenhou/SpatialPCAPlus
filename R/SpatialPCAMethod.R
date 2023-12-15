###################################################################################################
# Package: SpatialPCAplus
# Version: 0.0.9
# Date   : 2023-12-02
# Title : Spatial PCA with NNGP approximation
# Authors: Y. Hou
# Contacts: yirenhou@umich.edu
#          University of Michigan, Department of Biostatistics
###################################################################################################

#' SpatialPCAMethodEquations 
#'
#' @param kernelMatrix the kernel matrix K 
#' @param X scaled location matrix
#' @param tau parameter default is 1
#' @param equation test which equation
#'
#' @return equation output
#'
#' @export
SpatialPCAEquations_InverseKernel <- function(kernelMatrix, X, tau = 1, equation) {
  eigen_res = eigen(kernelMatrix)
  delta = eigen_res$values
  U = eigen_res$vectors
  n = ncol(U)
  
  eq = 0
  if (equation == 1) {
    #(M + tau^-1K^-1)^-1 
    part1 = U%*%solve((diag(1,n) + 1/tau * diag(1/delta)))%*%t(U)
    eq = part1 - part1%*%X%*%solve(-t(X)%*%X + t(X)%*%part1%*%X)%*%t(X)%*%part1
    
  } else if (equation == 2) {
    #|tau*K + In|
    eq = determinant(tau*diag(delta) + diag(1,n))
    
  } 
  return(eq)
}

#
#
#SpatialPCAMethod_InverseKernel - extracts inversion of kernel K in 
#SpatialPCA_EstimateLoading.R to compare with the inversion of 
#kernel K in NNGP_InverseKernel.R
#
#
SpatialPCA_InverseKernel <- function(kernelMatrix) {
  eigen_res = eigen(kernelMatrix)
  delta = eigen_res$values
  U = eigen_res$vectors
  
  inverseK = U%*%diag(1/delta)%*%t(U)
  return(inverseK)
}





