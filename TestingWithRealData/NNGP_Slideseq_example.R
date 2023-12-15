###################################################################################################
# Package: SpatialPCAplus
# Version: 0.0.9
# Date   : 2023-12-02
# Title : Spatial PCA with NNGP approximation
# Authors: Y. Hou
# Contacts: yirenhou@umich.edu
#          University of Michigan, Department of Biostatistics
###################################################################################################

library(ggplot2)
library(Rcpp)
library(SeuratObject)
library(Seurat)
library(ggpubr)
library(SpatialPCA)
source("nngp_approx.R")
source("SpatialPCA_EstimateLoading_NNGP.R")
source("SpatialPCA_SpatialPCs_NNGP.R")

load("slideseq.rds") 
print(dim(sp_count)) # The count matrix
print(dim(location)) # The location matrix



sp_count= sp_count[,1:10000]
location = as.matrix(location)[1:10000,]# location matrix: n x 2, count matrix: g x n.
# here n is spot number, g is gene number.
# here the column names of sp_count and rownames of location should be matched
slideseq = CreateSpatialPCAObject(counts=sp_count, location=location, project = "SpatialPCA",gene.type="spatial",sparkversion="sparkx",numCores_spark=1, customGenelist=NULL,min.loctions = 20, min.features=20)



slideseq.withKernel = SpatialPCA_buildKernel(slideseq, kerneltype="gaussian", bandwidthtype="Silverman",bandwidth.set.by.user=NULL,sparseKernel=FALSE,sparseKernel_tol=1e-20,sparseKernel_ncore=10)
slideseq = SpatialPCA_EstimateLoading(slideseq.withKernel,fast=TRUE,SpatialPCnum=20)
slideseq = SpatialPCA_SpatialPCs(slideseq, fast=TRUE)




system.time(SpatialPCA_EstimateLoading(slideseq.withKernel,fast=FALSE,SpatialPCnum=20))
system.time(SpatialPCA_SpatialPCs(slideseq, fast=FALSE))



# NNGP method
nngpOutput = nngp_approx(y = count_sub, m = 25, location = t(location),kernelMatrix = slideseq.withKernel@kernelmat)

slideseq.test = SpatialPCA_EstimateLoading_NNGP(slideseq.withKernel,fast = FALSE, SpatialPCnum = 20,
                                            sigma_inverse = nngpOutput$sigma.inverse, 
                                            sigma_inverse_det = nngpOutput$sigma.inverse.det, 
                                            sigma_det = nngpOutput$sigma.det)

slideseq.test = SpatialPCA_SpatialPCs_NNGP(slideseq.test, fast=FALSE,
                                       sigma_inverse = nngpOutput$sigma.inverse, 
                                       sigma_inverse_det = nngpOutput$sigma.inverse.det, 
                                       sigma_det = nngpOutput$sigma.det)

system.time(nngp_approx(y = rawcount, m = 25, location = t(location),kernelMatrix = slideseq.withKernel@kernelmat))
system.time(SpatialPCA_EstimateLoading_NNGP(slideseq.withKernel,fast = FALSE, SpatialPCnum = 20,
                                            sigma_inverse = nngpOutput$sigma.inverse, 
                                            sigma_inverse_det = nngpOutput$sigma.inverse.det, 
                                            sigma_det = nngpOutput$sigma.det))
system.time(SpatialPCA_SpatialPCs_NNGP(slideseq.test, fast=FALSE,
                                       sigma_inverse = nngpOutput$sigma.inverse, 
                                       sigma_inverse_det = nngpOutput$sigma.inverse.det, 
                                       sigma_det = nngpOutput$sigma.det))


   
slideseq = slideseq.test


clusterlabel= louvain_clustering(clusternum=8,latent_dat=slideseq@SpatialPCs,knearest=round(sqrt(dim(slideseq@SpatialPCs)[2])) )
cbp=c("#9C9EDE" ,"#5CB85C" ,"#E377C2", "#4DBBD5" ,"#FED439" ,"#FF9896", "#FFDC91")
plot_cluster(location=xy_coords,clusterlabel_refine,pointsize=1.5,text_size=20 ,title_in=paste0("SpatialPCA"),color_in=cbp,legend="right")

