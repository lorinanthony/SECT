### Clear Console ###
cat("\014")

### Clear Environment ### 
rm(list = ls(all = TRUE))

### Load in the R libraries ###
library(BGLR)
library(doParallel)
library(gdata)
library(Rcpp)
library(RcppArmadillo)
library(RcppParallel)
library(R.matlab)

######################################################################################
######################################################################################
######################################################################################

### Run the MATLAB Code that creates the EC Matrices for each image ###
#system("/Users/lorincrawford/Applications/MATLAB_R2017a.app/bin/matlab -nodisplay -r \"run('~/Dropbox/Columbia Radiogenomics/Software/CompEC.m'); exit\"")

### Load in the structural array that holds the EC Matrices ###
#Shapes = readMat('~/Dropbox/Columbia Radiogenomics/Data/ECs.mat')
#ECs = matrix(unlist(Shapes$Shapes[seq(2,length(Shapes$Shapes),2)]),nrow = length(Shapes$Shapes)/2,byrow = TRUE)

######################################################################################
######################################################################################
######################################################################################

### Run the R Code that creates the EC Matrices for each image ###
source("~/Dropbox/Columbia Radiogenomics/Software/EC3D.R")

### Set up the Parameters ###
setwd("~/Dropbox/Columbia Radiogenomics/Data")
startdir = getwd()
in.dir = "MITKSegmentations"
out.file = "~/Dropbox/Columbia Radiogenomics/Data/MRIECs.RData"

### Run The Euler Characteristic Function ###
ecf = ecf(in.dir = dir,out.file = out.file,img.dir = "baseline/Segmentations/enh",first.only = FALSE)

### Load in the List that holds the Euler Characteristic (EC) Curves for the TCIA Samples ###
load('~/Dropbox (Personal)/Columbia Radiogenomics/Data/MRI_ECs.RData')
nrot = ncol(MRI_list[[1]]$EC); stepsize = nrow(MRI_list[[1]]$EC)
ECs = matrix(nrow = length(MRI_list),ncol = nrot*stepsize)
rownames(ECs) = 1:nrow(ECs)
dim(ECs)

### Place the Curves in an nxp Matrix with Patient Tags as the Row Names ###
for(i in 1:nrow(ECs)){
  ECs[i,] = c(MRI_list[[i]]$EC)
  rownames(ECs)[i] = MRI_list[[i]]$name
}

### Remove the Vectors of Zeros where a New MRI Slice Begins ###
ECs = ECs[,-seq(101,ncol(ECs),by=101)]

