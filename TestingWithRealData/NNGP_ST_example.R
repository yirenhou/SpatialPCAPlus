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


load("Tumor_data.RData") 
print(dim(rawcount)) # The count matrix
print(dim(location)) # The location matrix


y = rawcount
m = 25
location = t(location)



# location matrix: n x 2, count matrix: g x n.
# here n is spot number, g is gene number.
# here the column names of sp_count and rownames of location should be matched
ST = CreateSpatialPCAObject(counts=rawcount, location=location, project = "SpatialPCA",gene.type="spatial",sparkversion="sparkx", gene.number=3000,customGenelist=NULL,min.loctions = 20, min.features=20)


ST.withKernel = SpatialPCA_buildKernel(ST, kerneltype="gaussian", bandwidthtype="SJ")
ST = SpatialPCA_EstimateLoading(ST.withKernel,fast=FALSE,SpatialPCnum=20)
ST = SpatialPCA_SpatialPCs(ST, fast=FALSE)

kernelMatrix = ST.withKernel@kernelmat
system.time(SpatialPCA_EstimateLoading(ST.withKernel,fast=FALSE,SpatialPCnum=20))
system.time(SpatialPCA_SpatialPCs(ST, fast=FALSE))



# NNGP method
nngpOutput = nngp_approx(y = rawcount, m = 25, location = t(location),kernelMatrix = ST.withKernel@kernelmat)

ST.test = SpatialPCA_EstimateLoading_NNGP(ST.withKernel,fast = FALSE, SpatialPCnum = 20,
                                          sigma_inverse = nngpOutput$sigma.inverse, 
                                          sigma_inverse_det = nngpOutput$sigma.inverse.det, 
                                          sigma_det = nngpOutput$sigma.det)

ST.test = SpatialPCA_SpatialPCs_NNGP(ST.test, fast=FALSE,
                                     sigma_inverse = nngpOutput$sigma.inverse, 
                                     sigma_inverse_det = nngpOutput$sigma.inverse.det, 
                                     sigma_det = nngpOutput$sigma.det)

system.time(nngp_approx(y = rawcount, m = 25, location = t(location),kernelMatrix = ST@kernelmat))
system.time(SpatialPCA_EstimateLoading_NNGP(ST.withKernel,fast = FALSE, SpatialPCnum = 20,
                                            sigma_inverse = nngpOutput$sigma.inverse, 
                                            sigma_inverse_det = nngpOutput$sigma.inverse.det, 
                                            sigma_det = nngpOutput$sigma.det))
system.time(SpatialPCA_SpatialPCs_NNGP(ST.test, fast=FALSE,
                                       sigma_inverse = nngpOutput$sigma.inverse, 
                                       sigma_inverse_det = nngpOutput$sigma.inverse.det, 
                                       sigma_det = nngpOutput$sigma.det))



###When using applied NNGP method###
ST = ST.test


##Detect spatial domain
clusterlabel= walktrap_clustering(7, ST.test@SpatialPCs,round(sqrt(dim(ST.test@location)[1])))
clusterlabel_refine=refine_cluster_10x(clusterlabel,ST.test@location,shape="square")

##Spatial domains detected by SpatialPCA
# set color
cbp_spatialpca = c(  "mediumaquamarine", "chocolate1","dodgerblue",  "#F0E442","palegreen4","lightblue2","plum1")
# visualize the cluster
plot_cluster(legend="right",location=ST.test@location,clusterlabel_refine,pointsize=5,text_size=20 ,title_in=paste0("SpatialPCA"),color_in=cbp_spatialpca)



##Trajectory inference
library(slingshot)
# trajectory on the whole tissue slice
sim = SingleCellExperiment(assays = rawcount)
reducedDims(sim) = SimpleList(DRM = t(ST@SpatialPCs))
colData(sim)$Walktrap = factor(clusterlabel_refine)    
# in this data we set tumor region as start cluster
sim  =slingshot(sim, clusterLabels = 'Walktrap', reducedDim = 'DRM',start.clus="2" ) 


# focus on tumor and its surrounding region
tumor_ind = which(clusterlabel_refine %in% c(2,3,7))
sim_tumor = SingleCellExperiment(assays = rawcount[,tumor_ind])
reducedDims(sim_tumor) = SimpleList(DRM = t(ST@SpatialPCs[,tumor_ind]))
colData(sim_tumor)$Walktrap = factor(clusterlabel_refine[tumor_ind])    
sim_tumor  =slingshot(sim_tumor, clusterLabels = 'Walktrap', reducedDim = 'DRM',start.clus="2" ) 
# in this data we set tumor region as start cluster
summary(sim_tumor@colData@listData)

# visualize on whole tissue
pseudotime_traj1_tumor = sim_tumor@colData@listData$slingPseudotime_1
clusterlabels_tumor = clusterlabel_refine[tumor_ind]
tumor_ind = which(clusterlabel_refine %in% c(2,3,7))
pseudotime_traj1 = sim@colData@listData$slingPseudotime_1
pseudotime_traj1[-tumor_ind]=NA
pseudotime_traj1[tumor_ind]=pseudotime_traj1_tumor

gridnum = 10
color_in=c(  "mediumaquamarine", "chocolate1","dodgerblue",  "#F0E442","palegreen4","lightblue2","plum1","black","#CC79A7","mediumpurple","seagreen1")

p_traj1 = plot_trajectory(pseudotime_traj1, location,clusterlabel,gridnum,color_in,pointsize=5 ,arrowlength=0.3,arrowsize=1.3,textsize=15 )
print(ggarrange( p_traj1[[4]],p_traj1[[1]],
                 ncol = 2, nrow = 1))


##High resolution map reconstruction
STsimu_high_ST = SpatialPCA_highresolution(ST, platform="ST",newlocation=NULL)
cluster_SpatialPCA_high = walktrap_clustering(7, latent_dat=STsimu_high_ST@highPCs,200)
color_in=c(  "palegreen4", "chocolate1","plum1",  "#F0E442","mediumaquamarine","dodgerblue","lightblue2")
title_in="SpatialPCA High resolution"
plot_cluster(STsimu_high_ST@highPos, as.character(cluster_SpatialPCA_high), pointsize=2,text_size=20 ,title_in,color_in,legend="bottom")


