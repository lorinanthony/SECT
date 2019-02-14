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

######################################################################################
######################################################################################
######################################################################################

### Specify a Function for Standard Error ###
error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
  if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
    stop("vectors must be same length")
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, lwd = 2, ...)
}

######################################################################################
######################################################################################
######################################################################################

### Change the Working Directory ###
setwd("~/Dropbox/Columbia Radiogenomics/Analysis/Results/")

### Load in the Cross-Validation Results ###
load("Linear_SECT.RData"); A = Res
load("Gauss_SECT.RData"); B = Res
load("Cauchy_SECT.RData"); C = Res

######################################################################################
######################################################################################
######################################################################################

### Analysis of Methods and Data According to Squared Correlation Coefficient (R^2) ###

### Set the Number of Splits Used in the Simulation ###
n.splits = nrow(A[[1]]); FMR = C

### Call the Phenotypes ###
DFS = FMR[[2]][,5:8]
OS = FMR[[1]][,5:8]

### Disease Free Survival ###
countCOR = matrix(0,nrow = n.splits,ncol = 4)
for(i in 1:nrow(FMR[[1]])){
  countCOR[i,which(DFS[i,] == max(DFS[i,]))] = 1
}

colMeans(countCOR)

### Overall Survival ###
countCOR = matrix(0,nrow = n.splits,ncol = 4)
for(i in 1:nrow(FMR[[1]])){
  countCOR[i,which(OS[i,] == max(OS[i,]))] = 1
}

colMeans(countCOR)

### Get the Overall Mean and Standard Error of the Results ###
round(cbind(colMeans(DFS),apply(DFS,2,sd)/sqrt(n.splits)),3)
round(cbind(colMeans(OS),apply(OS,2,sd)/sqrt(n.splits)),3)

######################################################################################
######################################################################################
######################################################################################

### Set the Working Directory for Figures ###
setwd("~/Dropbox/Columbia Radiogenomics/Analysis/Figures")

### Look at the Squared Correlation Coefficient ###
dfs.cor = c(colMeans(A[[2]][,5:8]),colMeans(B[[2]][,5:8]),colMeans(C[[2]][,5:8]))
os.cor = c(colMeans(A[[1]][,5:8]),colMeans(B[[1]][,5:8]),colMeans(C[[1]][,5:8]))

dfs.cor.sd = c(apply(A[[2]][,5:8],2,sd),apply(B[[2]][,5:8],2,sd),apply(C[[2]][,5:8],2,sd))/10
os.cor.sd = c(apply(A[[1]][,5:8],2,sd),apply(B[[1]][,5:8],2,sd),apply(C[[1]][,5:8],2,sd))/10
