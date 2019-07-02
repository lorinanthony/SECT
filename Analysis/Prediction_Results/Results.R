### Load in the R libraries ###
library(BGLR)
library(doParallel)
library(gdata)
library(Rcpp)
library(RcppArmadillo)
library(RcppParallel)

######################################################################################
######################################################################################
######################################################################################

### Load in the Prediction Results ###
load("./Analysis/Prediction_Results/Linear_SECT.RData"); A = Res
load("./Analysis/Prediction_Results/Gauss_SECT.RData"); B = Res
load("./Analysis/Prediction_Results/Cauchy_SECT.RData"); C = Res

### Take the Number of Testing Splits ###
n.splits = nrow(A[[1]])

######################################################################################
######################################################################################
######################################################################################

### Find the Optimal Method Count: Linear Covariance Function ###
FMR = A

### Call the Phenotypes ###
DFS = FMR[[2]]
OS = FMR[[1]]

### Disease Free Survival ###
A_opt_DFS = matrix(0,nrow = n.splits,ncol = 4)
for(i in 1:nrow(FMR[[1]])){
  A_opt_DFS[i,which(DFS[i,] == max(DFS[i,]))] = 1
}

### Overall Survival ###
A_opt_OS = matrix(0,nrow = n.splits,ncol = 4)
for(i in 1:nrow(FMR[[1]])){
  A_opt_OS[i,which(OS[i,] == max(OS[i,]))] = 1
}

######################################################################################
######################################################################################
######################################################################################

### Find the Optimal Method Count: Gaussian Covariance Function ###
FMR = B

### Call the Phenotypes ###
DFS = FMR[[2]]
OS = FMR[[1]]

### Disease Free Survival ###
B_opt_DFS = matrix(0,nrow = n.splits,ncol = 4)
for(i in 1:nrow(FMR[[1]])){
  B_opt_DFS[i,which(DFS[i,] == max(DFS[i,]))] = 1
}

### Overall Survival ###
B_opt_OS = matrix(0,nrow = n.splits,ncol = 4)
for(i in 1:nrow(FMR[[1]])){
  B_opt_OS[i,which(OS[i,] == max(OS[i,]))] = 1
}

######################################################################################
######################################################################################
######################################################################################

### Find the Optimal Method Count: Cauchy Covariance Function ###
FMR = C

### Call the Phenotypes ###
DFS = FMR[[2]]
OS = FMR[[1]]

### Disease Free Survival ###
C_opt_DFS = matrix(0,nrow = n.splits,ncol = 4)
for(i in 1:nrow(FMR[[1]])){
  C_opt_DFS[i,which(DFS[i,] == max(DFS[i,]))] = 1
}

### Overall Survival ###
C_opt_OS = matrix(0,nrow = n.splits,ncol = 4)
for(i in 1:nrow(FMR[[1]])){
  C_opt_OS[i,which(OS[i,] == max(OS[i,]))] = 1
}

######################################################################################
######################################################################################
######################################################################################

### Load in the Cross-Validated Bandwidths ###
theta = seq(from = 0.1, to = 10, by = 0.1)

### Call from the Spectrum of Considered Bandwidths ###
load("./Analysis/Cross_Validation_Results/GaussCV_Results.RData");
gauss.theta.dfs = c(theta[cvs[[2]][1]],theta[cvs[[2]][3]],theta[cvs[[2]][2]],theta[cvs[[2]][4]])
gauss.theta.os = c(theta[cvs[[1]][1]],theta[cvs[[1]][3]],theta[cvs[[1]][2]],theta[cvs[[1]][4]])

### Call from the Spectrum of Considered Bandwidths ###
load("./Analysis/Cross_Validation_Results/CauchyCV_Results.RData");
cauchy.theta.dfs = c(theta[cvs[[2]][1]],theta[cvs[[2]][3]],theta[cvs[[2]][2]],theta[cvs[[2]][4]])
cauchy.theta.os = c(theta[cvs[[1]][1]],theta[cvs[[1]][3]],theta[cvs[[1]][2]],theta[cvs[[1]][4]])

######################################################################################
######################################################################################
######################################################################################

### Present Results for Linear Covariance Function ###

### Disease Free Survival (DFS) ###
mat_DFS = matrix(NA,nrow = 4,ncol = 4)
colnames(mat_DFS) = c("Mean R^2 ","SD R^2","Optimal%","Bandwidth")
rownames(mat_DFS) = c("Gene Expression","Morphometrics","Geometrics","SECT")

mat_DFS[,1] = round(colMeans(A[[2]],na.rm = TRUE),3); mat_DFS[,2] = round(apply(A[[2]],2,sd,na.rm = TRUE),2)/10
mat_DFS[,3] = colMeans(A_opt_DFS)

### Overall Survival ###
mat_OS = matrix(NA,nrow = 4,ncol = 4)
colnames(mat_OS) = c("Mean R^2 ","SD R^2","Optimal%","Bandwidth")
rownames(mat_OS) = c("Gene Expression","Morphometrics","Geometrics","SECT")

mat_OS[,1] = round(colMeans(A[[1]],na.rm = TRUE),3); mat_OS[,2] = round(apply(A[[1]],2,sd,na.rm = TRUE),2)/10
mat_OS[,3] = colMeans(A_opt_OS)

### Present the Results ###
cat("\nPredictive Results using the Linear Covariance Function.\n")
cat("\nDisease Free Survival (DFS):\n");print(mat_DFS)
cat("\nOverall Survival (OS):\n"); print(mat_OS)

######################################################################################
######################################################################################
######################################################################################

### Present Results for Gaussian Covariance Function ###

### Disease Free Survival (DFS) ###
mat_DFS = matrix(NA,nrow = 4,ncol = 4)
colnames(mat_DFS) = c("Mean R^2 ","SD R^2","Optimal%","Bandwidth")
rownames(mat_DFS) = c("Gene Expression","Morphometrics","Geometrics","SECT")

mat_DFS[,1] = round(colMeans(B[[2]],na.rm = TRUE),3); mat_DFS[,2] = round(apply(B[[2]],2,sd,na.rm = TRUE),2)/10
mat_DFS[,3] = colMeans(B_opt_DFS); mat_DFS[,4] = gauss.theta.dfs

### Overall Survival ###
mat_OS = matrix(NA,nrow = 4,ncol = 4)
colnames(mat_OS) = c("Mean R^2 ","SD R^2","Optimal%","Bandwidth")
rownames(mat_OS) = c("Gene Expression","Morphometrics","Geometrics","SECT")

mat_OS[,1] = round(colMeans(B[[1]],na.rm = TRUE),3); mat_OS[,2] = round(apply(B[[1]],2,sd,na.rm = TRUE),2)/10
mat_OS[,3] = colMeans(B_opt_OS); mat_OS[,4] = gauss.theta.os

### Present the Results ###
cat("\nPredictive Results using the Gaussian Covariance Function.\n")
cat("\nDisease Free Survival (DFS):\n");print(mat_DFS)
cat("\nOverall Survival (OS):\n");print(mat_OS)

######################################################################################
######################################################################################
######################################################################################

### Present Results for Cauchy Covariance Function ###

### Disease Free Survival (DFS) ###
mat_DFS = matrix(NA,nrow = 4,ncol = 4)
colnames(mat_DFS) = c("Mean R^2 ","SD R^2","Optimal%","Bandwidth")
rownames(mat_DFS) = c("Gene Expression","Morphometrics","Geometrics","SECT")

mat_DFS[,1] = round(colMeans(C[[2]],na.rm = TRUE),3); mat_DFS[,2] = round(apply(C[[2]],2,sd,na.rm = TRUE),2)/10
mat_DFS[,3] = colMeans(C_opt_DFS); mat_DFS[,4] = cauchy.theta.dfs

### Overall Survival ###
mat_OS = matrix(NA,nrow = 4,ncol = 4)
colnames(mat_OS) = c("Mean R^2 ","SD R^2","Optimal%","Bandwidth")
rownames(mat_OS) = c("Gene Expression","Morphometrics","Geometrics","SECT")

mat_OS[,1] = round(colMeans(C[[1]],na.rm = TRUE),3); mat_OS[,2] = round(apply(C[[1]],2,sd,na.rm = TRUE),2)/10
mat_OS[,3] = colMeans(C_opt_OS); mat_OS[,4] = cauchy.theta.os

### Present the Results ###
cat("\nPredictive Results using the Cauchy Covariance Function.\n")
cat("\nDisease Free Survival (DFS):\n");print(mat_DFS)
cat("\nOverall Survival (OS):\n");print(mat_OS)