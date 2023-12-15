###############################################################################################################
# Package: SpatialPCAplus
# Version: 0.0.9
# Date   : 2023-12-02
# Title : Spatial PCA with NNGP approximation
# Authors: Y. Hou
# Contacts: yirenhou@umich.edu
#          University of Michigan, Department of Biostatistics
###########################################################################################################

#' @title nngp_approx
#' @name nngp_approx
#'
#' @param kernelMatrix the kernel matrix K from SpatialPCA
#' @param y gene expression matrix
#' @param m number of nearest number
#' @param location A n by d matrix of cell/spot location coordinates.
#
#' @return sigma.inverse, sigma.inverse.det, sigma.det, F.vector
#' @import foreach 
#' @import Matrix
#'
#' @export
library(foreach)
library(Matrix)
nngp_approx <- function(y, m, location, kernelMatrix) {
  if(is.vector(y)) {
    n = length(y)
  } else {
    n = ncol(y)
  }
  
  k = nrow(location)
  x.coord = location[1,]
  ##working with directed acyclic graphs (DAGS) to define each set of neighbors
  #order input points si's according to some criterion (sort by first coordinate)
  ordered.loc = t(location[,order(x.coord)])
  
  #define a valid distance metric between points (Euclidean distance)
  euc.dist <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))
  dist.si = outer(1:n,1:n, Vectorize(function(i, j) sqrt(sum((ordered.loc[i,] - ordered.loc[j,])^2))))
  
  #for each element i, define m nearest neighbors N(si) to be the set of points sj1,...,sjm
  #such that all neighbor indices j1,...,jm are less than i
  B.matrix = c()
  F.vector = c()
  i = 1
  for (i in 1:n) {
    which.si = which(location[1,] == ordered.loc[i,1] & location[2,] == ordered.loc[i,2])
    if ((i - m) == (1 - m)) {
      Ksi = kernelMatrix[which.si,which.si]
      
      ##densities for NNGP
      #Bsi = Ksi,N(si) %*% KN(si)
      Bsi = 0
      
      #Fsi = Ksi - Ksi,N(si) %*% solve(KN(si)) %*% t(Ksi,N(si))
      Fsi = Ksi
      
      ##contribution of each sample xj to the term Bxi
      bstar = rep(0,n)
      bstar[i] = 1
      
      #construct Bsi* and B
      Bstar = as.matrix(bstar)
      
      B.matrix = cbind(B.matrix,Bstar)
      F.vector = c(F.vector,Fsi)
    } else if ((i - m) == (2 - m)) {
      neighbors.si = ordered.loc[1,]
      which.Nsi = which(location[1,] == neighbors.si[1] & location[2,] == neighbors.si[2])
      
      ##construct Kernel Matrix using elements in N(si)
      Ksi = kernelMatrix[which.si,which.si]
      KNsi = kernelMatrix[which.Nsi,which.Nsi]
      Ksi.Nsi = kernelMatrix[which.si,which.Nsi]
      
      ##densities for NNGP
      #Bsi = Ksi,N(si) %*% KN(si)
      Bsi = Ksi.Nsi %*% KNsi
      #Fsi = Ksi - Ksi,N(si) %*% solve(KN(si)) %*% t(Ksi,N(si))
      Fsi = Ksi - Ksi.Nsi%*%solve(KNsi)%*%Ksi.Nsi
      
      
      ##contribution of each sample xj to the term Bxi
      bstar = rep(0,n)
      bstar[i] = 1
      bstar[1] = Bsi
      
      #construct Bsi* and B
      Bstar = as.matrix(bstar)
      
      B.matrix = cbind(B.matrix,Bstar)
      F.vector = c(F.vector,Fsi)
    } else if ((i - m) < 1 & (i - m) > (2 - m)) {
      neighbors.si = ordered.loc[c(1:(i-1)),]
      
      which.Nsi = sapply(1:nrow(neighbors.si), function(j) which(location[1,] == neighbors.si[j,1] & location[2,] == neighbors.si[j,2]))
      
      ##yN(si) is the subvector of y containing all yj such that sj is in N(si)
      #yNsi = neighbors.si[,k+1]
      
      ##construct Kernel Matrix using elements in N(si)
      Ksi = kernelMatrix[which.si,which.si]
      KNsi = kernelMatrix[which.Nsi,which.Nsi]
      Ksi.Nsi = kernelMatrix[which.si,which.Nsi]
      
      ##densities for NNGP
      #Bsi = Ksi,N(si) %*% KN(si)
      Bsi = Ksi.Nsi %*% KNsi
      #Fsi = Ksi - Ksi,N(si) %*% solve(KN(si)) %*% t(Ksi,N(si))
      Fsi = Ksi - Ksi.Nsi%*%solve(KNsi)%*%Ksi.Nsi
      
      
      ##contribution of each sample xj to the term Bxi
      bstar = rep(0,n)
      bstar[i] = 1
      bstar[c(1:(i-1))] = Bsi
      
      #construct Bsi* and B
      Bstar = as.matrix(bstar)
      
      B.matrix = cbind(B.matrix,Bstar)
      F.vector = c(F.vector,Fsi)
    } else {
      neighbor.order = order(dist.si[i,1:(i - 1)])[1:m]
      neighbors.si = ordered.loc[sort(neighbor.order),]
      
      which.Nsi = sapply(1:nrow(neighbors.si), function(j) which(location[1,] == neighbors.si[j,1] & location[2,] == neighbors.si[j,2]))
      
      ##yN(si) is the subvector of y containing all yj such that sj is in N(si)
      #yNsi = neighbors.si[,k+1]
      
      ##construct Kernel Matrix using elements in N(si)
      Ksi = kernelMatrix[which.si,which.si]
      KNsi = kernelMatrix[which.Nsi,which.Nsi]
      Ksi.Nsi = kernelMatrix[which.si,which.Nsi]
      
      ##densities for NNGP
      #Bsi = Ksi,N(si) %*% KN(si)
      Bsi = Ksi.Nsi %*% KNsi
      
      #Fsi = Ksi - Ksi,N(si) %*% solve(KN(si)) %*% t(Ksi,N(si))
      Fsi = Ksi - Ksi.Nsi%*%solve(KNsi)%*%Ksi.Nsi
      
      
      ##contribution of each sample xj to the term Bxi
      bstar = rep(0,n)
      bstar[i] = 1
      bstar[neighbor.order] = Bsi
      
      #construct Bsi* and B
      Bstar = as.matrix(bstar)
      
      B.matrix = cbind(B.matrix,Bstar)
      F.vector = c(F.vector,Fsi)
    }
    
  }
  
  B.matrix = t(B.matrix)
  
  ##define F, nxn diagonal matrix
  F.matrix.inverse = diag(1 / F.vector)
  
  ##precision matrix (inverse of covariance matrix) given by t(B) %*% solve(F) %*% B
  sigma.inverse = t(B.matrix)%*%F.matrix.inverse%*%B.matrix
  
  ##return precision matrix
  return(list(sigma.inverse = sigma.inverse, sigma.inverse.det = prod(1/F.vector), sigma.det = prod(F.vector), F.vector = F.vector))
}


