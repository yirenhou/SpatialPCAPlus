# SpatialPCAPlus
This contains all the codes for BIOSTAT 615 Project SpatialPCAPlus. Below describes each folder and code for easier navigation.
#
Data folder - contains cortex data from DLPFC (LIBD), Mouse cerebellum data by Slide-seq, HER2 tumpor data used for real data analysis.
#
NewTestingCodes folder - responsible for test cases comparing computational time between NNGP approximation of inverse matrix K and eigen-decomposition for inverse
                  matrix in SpatialPCA. Utilized functions BuildTestingKernels (BuildTestingKernels.R), nngp_approx (NNGP_InverseKernel.R), nngpMethodEquations_InverseKernel (NNGPMethod.R), SpatialPCAEquations_InverseKernel and SpatialPCA_InverseKernel (SpatialPCAMethod.R)
#
R folder - contains all the R files with functions other testing and analysis files need.
#
TestingWithRealData folder - responsible for real data analysis. Warning: DLPFC and Slide-seq are very large dataset which require long computation time. Utilized functions nngp_approx, SpatialPCA_EstimateLoading_NNGP, and SpatialPCA_SpatialPCs_NNGP.
