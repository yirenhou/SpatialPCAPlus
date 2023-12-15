###################################################################################################
# Package: SpatialPCAplus
# Version: 0.0.9
# Date   : 2023-12-02
# Title : Spatial PCA with NNGP approximation
# Authors: Y. Hou
# Contacts: yirenhou@umich.edu
#          University of Michigan, Department of Biostatistics
###################################################################################################

#' NNGP equations
#'
#'
#' @param kernelMatrix the kernel matrix K 
#' @param X scaled location matrix
#' @param F vector from Nearest Neighbor Gaussian Process 
#' @param tau parameter default is 1
#' @param equation test which equation
#'
#' @return equation output
#'
#' @export
nngpMethodEquations_InverseKernel <- function(inverseMatrix, F.vector, X, tau = 1,equation) {
  n = ncol(inverseMatrix)
  
  eq = 0
  if (equation == 1) {
    #(M + tau^-1K^-1)^-1 
    eq = solve(diag(1,n) + 1/tau*inverseMatrix + X%*%(-solve(t(X)%*%X))%*%t(X))

  } else if (equation == 2) {
    #|tau*K + In|
    eq = prod((tau*F.vector + 1))
      
  } 
    
  return(eq)
}
