### Clear Console ###
cat("\014")

### Clear Environment ### 
rm(list = ls(all = TRUE))

### Load in the R libraries ###
library(BGLR)
library(doParallel)
library(Rcpp)
library(RcppArmadillo)
library(RcppParallel)
#library(R.matlab)

######################################################################################
######################################################################################
######################################################################################

### Run the MATLAB Code that creates the EC Matrices for each image ###
#system("~/Applications/MATLAB_R2017a.app/bin/matlab -nodisplay -r \"run('./Software/CompEC.m'); exit\"")

### Load in the structural array that holds the EC Matrices ###
#Shapes = readMat('./Data/ECs.mat')
#ECs = matrix(unlist(Shapes$Shapes[seq(2,length(Shapes$Shapes),2)]),nrow = length(Shapes$Shapes)/2,byrow = TRUE)

######################################################################################
######################################################################################
######################################################################################

### Run the R Code that creates the EC Matrices for each image ###
source("./Software/EC3D.R")

ptm <- proc.time() #Start clock for timing

### Set up the Parameters ###
startdir = "./Data"
in.dir = "./Data/MITKSegmentations"
out.file = "./Data/MRIECs.RData"

### Run The Euler Characteristic Function ###
ecf = ecf(in.dir = in.dir,out.file = out.file,img.dir = "baseline/Segmentations/enh",first.only = FALSE)

### Load in the List that holds the Euler Characteristic (EC) Curves for the TCIA Samples ###
load("./Data/MRI_ECs.RData")
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

proc.time() - ptm #Stop clock and check

######################################################################################
######################################################################################
######################################################################################

### Plot the SECT for Different Samples in 1 Direction ###
plot(ECs[1,1:100],type = "l", col = "blue", lwd = 2,bty = "n",ylab = "Smooth Euler Characteristic Transform (SECT)",xlab = "Sublevel Sets", ylim = c(1.5,4))
lines(ECs[2,1:100],col = "red", lty = 2, lwd = 2)
lines(ECs[4,1:100],col = "forest green", lty = 4, lwd = 2)
legend("top",legend = c("Patient #1","Patient #2","Patient #3"),col = c("blue","red","forest green"),lty = c(1,2,4),lwd = 2,horiz = TRUE,bty = "n")
