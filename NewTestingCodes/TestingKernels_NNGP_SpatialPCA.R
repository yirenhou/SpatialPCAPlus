###################################################################################################
# Package: SpatialPCAplus
# Version: 0.0.9
# Date   : 2023-12-02
# Title : Spatial PCA with NNGP approximation
# Authors: Y. Hou
# Contacts: yirenhou@umich.edu
#          University of Michigan, Department of Biostatistics
###################################################################################################
#
#Testing different type of gaussian kernels for kernel inversion using NNGP or eigen-decom
#position used in SpatialPCA
#
#
#
source("BuildTestingKernels.R")
source("nngp_approx.R")
source("NNGPMethod.R")
source("SpatialPCAMethod.R")

###Random generated data from Uniform distribution
##5000 x 5000 kernel matrix test case
n = 5000
x.coord = runif(n,0,20)
y.coord = rnorm(n,0,25)
y = matrix(sample(2,n*n,replace = TRUE), ncol = n, nrow = n)

location = cbind(x.coord,y.coord)
X = scale(location)
kernelTest1 = BuildTestingKernels(location,bandwidth = 1)


system.time(nngp_approx(y = y, m = 30, location = t(location),kernelMatrix = kernelTest1))
system.time(SpatialPCA_InverseKernel(kernelTest1))


nngp_output = nngp_approx(y = y, m = 30, location = t(location),kernelMatrix = kernelTest1)
sPCA_output = SpatialPCA_InverseKernel(kernelTest1)


system.time(nngpMethodEquations_InverseKernel(nngp_output$sigma.inverse,nngp_output$F.vector,X = X,equation = 1))
system.time(SpatialPCAEquations_InverseKernel(kernelTest1,X = X, equation = 1))

##7000 x7000 kernel matrix test case
n = 7000
x.coord = runif(n,0,20)
y.coord = rnorm(n,0,25)
y = matrix(sample(2,n*n,replace = TRUE), ncol = n, nrow = n)

location = cbind(x.coord,y.coord)
X = scale(location)
kernelTest2 = BuildTestingKernels(location,bandwidth = 1)

system.time(nngp_approx(y = y, m = 30, location = t(location),kernelMatrix = kernelTest2))
system.time(SpatialPCA_InverseKernel(kernelTest2))

nngp_output = nngp_approx(y = y, m = 30, location = t(location),kernelMatrix = kernelTest2)
sPCA_output = SpatialPCA_InverseKernel(kernelTest2)


system.time(nngpMethodEquations_InverseKernel(nngp_output$sigma.inverse,nngp_output$F.vector,X = X,equation = 1))
system.time(SpatialPCAEquations_InverseKernel(kernelTest2,X = X, equation = 1))


###ST dataset
load("Tumor_data.RData") 
X = scale(location)
X = scale(location)

location = location[,1:2]
kernelTest3 = BuildTestingKernels(location,bandwidth = 1)

nngp_output = nngp_approx(y = rawcount, m = 30, location = t(location),kernelMatrix = kernelTest3)
sPCA_output = SpatialPCA_InverseKernel(kernelTest3)

system.time(nngp_approx(y = rawcount, m = 30, location = t(location),kernelMatrix = kernelTest3))
system.time(SpatialPCA_InverseKernel(kernelTest3))


nngp_output = nngp_approx(y = y, m = 30, location = t(location),kernelMatrix = kernelTest3)
sPCA_output = SpatialPCA_InverseKernel(kernelTest3)

system.time(nngpMethodEquations_InverseKernel(nngp_output$sigma.inverse,nngp_output$F.vector,X = X,equation = 1))
system.time(SpatialPCAEquations_InverseKernel(kernelTest3,X = X, equation = 1))

