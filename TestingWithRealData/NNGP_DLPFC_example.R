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


sample_names=c("151507", "151508", "151509", "151510", "151669", "151670", "151671" ,"151672","151673", "151674" ,"151675" ,"151676")
i=9 # Here we take the ith sample as example, in total there are 12 samples (numbered as 1-12), the user can test on other samples if needed.
clusterNum=c(7,7,7,7,5,5,5,5,7,7,7,7) # each sample has different ground truth cluster number
load(paste0("LIBD_sample",i,".RData")) 
print(dim(count_sub)) # The count matrix
print(dim(xy_coords)) # The x and y coordinates. We flipped the y axis for visualization.



# location matrix: n x 2, count matrix: g x n.
# here n is spot number, g is gene number.
count_sub = count_sub[1:20000,1:1500]
xy_coords = as.matrix(xy_coords)[1:1500,]
rownames(xy_coords) = colnames(count_sub) # the rownames of location should match with the colnames of count matrix
LIBD = CreateSpatialPCAObject(counts=count_sub, location=xy_coords, project = "SpatialPCA",gene.type="spatial",sparkversion="sparkx",numCores_spark=1,gene.number=3000, customGenelist=NULL,min.loctions = 20, min.features=20)



LIBD.withKernel = SpatialPCA_buildKernel(LIBD, kerneltype="gaussian", bandwidthtype="SJ",bandwidth.set.by.user=NULL)
LIBD = SpatialPCA_EstimateLoading(LIBD.withKernel,fast=FALSE,SpatialPCnum=20) 
LIBD = SpatialPCA_SpatialPCs(LIBD, fast=FALSE)




system.time(SpatialPCA_EstimateLoading(LIBD.withKernel,fast=FALSE,SpatialPCnum=20))
system.time(SpatialPCA_SpatialPCs(LIBD, fast=FALSE))



# NNGP method
nngpOutput = nngp_approx(y = count_sub, m = 25, location = t(xy_coords),kernelMatrix = LIBD@kernelmat)

LIBD.test = SpatialPCA_EstimateLoading_NNGP(LIBD.withKernel,fast = FALSE, SpatialPCnum = 20,
                                          sigma_inverse = nngpOutput$sigma.inverse, 
                                          sigma_inverse_det = nngpOutput$sigma.inverse.det, 
                                          sigma_det = nngpOutput$sigma.det)

LIBD.test = SpatialPCA_SpatialPCs_NNGP(LIBD.test, fast=FALSE,
                                     sigma_inverse = nngpOutput$sigma.inverse, 
                                     sigma_inverse_det = nngpOutput$sigma.inverse.det, 
                                     sigma_det = nngpOutput$sigma.det)

system.time(nngp_approx(y = count_sub, m = 25, location = t(xy_coords),kernelMatrix = LIBD.withKernel@kernelmat))
system.time(SpatialPCA_EstimateLoading_NNGP(LIBD.withKernel,fast = FALSE, SpatialPCnum = 20,
                                            sigma_inverse = nngpOutput$sigma.inverse, 
                                            sigma_inverse_det = nngpOutput$sigma.inverse.det, 
                                            sigma_det = nngpOutput$sigma.det))
system.time(SpatialPCA_SpatialPCs_NNGP(LIBD.test, fast=FALSE,
                                       sigma_inverse = nngpOutput$sigma.inverse, 
                                       sigma_inverse_det = nngpOutput$sigma.inverse.det, 
                                       sigma_det = nngpOutput$sigma.det))



LIBD = LIBD.test



clusterlabel= walktrap_clustering(clusternum=clusterNum[i],latent_dat=LIBD.test@SpatialPCs,knearest=70 ) 
# here for all 12 samples in LIBD, we set the same k nearest number in walktrap_clustering to be 70. 
# for other Visium or ST data, the user can also set k nearest number as round(sqrt(dim(SpatialPCAobject@SpatialPCs)[2])) by default.
clusterlabel_refine = refine_cluster_10x(clusterlabels=clusterlabel,location=LIBD.test@location,shape="hexagon")


cbp=c("#9C9EDE" ,"#5CB85C" ,"#E377C2", "#4DBBD5" ,"#FED439" ,"#FF9896", "#FFDC91")
plot_cluster(location=xy_coords,clusterlabel=clusterlabel_refine,pointsize=1.5,text_size=20 ,title_in=paste0("SpatialPCA"),color_in=cbp)



truth = KRM_manual_layers_sub$layer_guess_reordered[match(colnames(LIBD.test@normalized_expr),colnames(count_sub))]
cbp=c("#5CB85C" ,"#9C9EDE" ,"#FFDC91", "#4DBBD5" ,"#FF9896" ,"#FED439", "#E377C2", "#FED439")
plot_cluster(location=xy_coords,truth,pointsize=1.5,text_size=20 ,title_in=paste0("Ground truth"),color_in=cbp)



library(slingshot)
sim<- SingleCellExperiment(assays = count_sub)
reducedDims(sim) <- SimpleList(DRM = t(LIBD@SpatialPCs))
colData(sim)$clusterlabel <- factor(clusterlabel_refine)    
sim  <-slingshot(sim, clusterLabels = 'clusterlabel', reducedDim = 'DRM',start.clus="3" ) 
# in this data we set white matter region as start cluster, one can change to their preferred start region 

summary(sim@colData@listData)
pseudotime_traj1 = sim@colData@listData$slingPseudotime_1 # in this data only one trajectory was inferred
gridnum = 10
color_in = c("#9C9EDE" ,"#5CB85C" ,"#E377C2", "#4DBBD5" ,"#FED439" ,"#FF9896", "#FFDC91","black")
p_traj1 = plot_trajectory(pseudotime_traj1, LIBD@location,clusterlabel_refine,gridnum,color_in,pointsize=1.5 ,arrowlength=0.2,arrowsize=1,textsize=15 )
p_traj1$Arrowoverlay1
p_traj1$Pseudotime


