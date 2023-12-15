########################################################################################################################
# Package: SpatialPCA
# Version: 1.1.0
# Date   : 2021-10-27
# Title : Spatially Aware Dimension Reduction for Spatial Transcriptomics 
# Authors: L. Shang and X. Zhou
# Contacts: shanglu@umich.edu 
#          University of Michigan, Department of Biostatistics
####################################################################################################


#' Calculating Spatial PCs (latent factor matrix Z).
#' @param object SpatialPCA object.
#' @param fast Select fast=TRUE if the user wants to use low-rank approximation on the kernel matrix to calculate the spatial PCs, otherwise select FALSE. 
#' @param eigenvecnum When fast=TRUE, eigenvecnum is the number of top eigenvectors and eigenvalues to be used in low-rank approximation in the eigen decomposition step for kernel matrix. 
#' The default is NULL, if specified, it is recommended that these top eigen values explain >=90% of the variance. 
#' In estimating spatial PCs, we need larger number of eigenvectors in kernel matrix for more accurate estimation.
#' @return Returns SpatialPCA object with estimated Spatial PCs.
#' 
#' @import RSpectra
#' 
#' @export
SpatialPCA_SpatialPCs_NNGP = function(object,fast=FALSE,eigenvecnum=NULL,sigma_inverse, sigma_inverse_det, sigma_det){
  
  # suppressMessages(require(RSpectra))
  
  n = object@params$n
  PCnum = object@SpatialPCnum
  Z_hat = matrix(0, PCnum, n)
  tau = object@tau
  W_hat = object@W
  
  
  W_hat_t = t(W_hat)
  WtYM = W_hat_t%*% object@params$YM

  
  object@SpatialPCs = WtYM %*% solve(object@params$M + 1/tau*sigma_inverse)
  object@SpatialPCs = as.matrix(object@SpatialPCs)
  rm(W_hat_t)
  rm(WtYM)
  gc()
  
  
  return(object)
}

