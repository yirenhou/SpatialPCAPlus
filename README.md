# SpatialPCAPlus
This contains all the codes for BIOSTAT 615 Project SpatialPCAPlus. Below describes each folder and code for easier navigation.
#
Data folder - contains cortex data from DLPFC (LIBD), Mouse cerebellum data by Slide-seq, HER2 tumpor data used for real data analysis.
#
NewTestingCodes - responsible for test cases comparing computational time between NNGP approximation of inverse matrix K and eigen-decomposition for inverse
                  matrix in SpatialPCA
                - utilized functions BuildTestingKernels (BuildTestingKernels.R), nngp_approx (NNGP_InverseKernel.R), nngpMethodEquations_InverseKernel (NNGPMethod.R), SpatialPCAEquations_InverseKernel and SpatialPCA_InverseKernel (SpatialPCAMethod.R)

#

