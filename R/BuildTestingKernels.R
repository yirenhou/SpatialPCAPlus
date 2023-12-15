###################################################################################################
# Package: SpatialPCAplus
# Version: 0.0.9
# Date   : 2023-12-02
# Title : Spatial PCA with NNGP approximation
# Authors: Y. Hou
# Contacts: yirenhou@umich.edu
#          University of Michigan, Department of Biostatistics
###################################################################################################

#' 
#' Function builds Gaussian kernels for testing.
#' 
#' 
#' @param bandwidth A numeric value of bandwidth
#' @param location A n by d matrix of cell/spot location coordinates.
#' @return Gaussian kernel
#' @export
BuildTestingKernels <- function(location, bandwidth) {
  K = exp(-1*as.matrix(dist(location)^2)/bandwidth)
  return(K)
}


