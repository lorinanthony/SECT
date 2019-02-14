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

### Load in the C++ BAKR Kernel functions ###
sourceCpp("~/Dropbox (Personal)/Columbia Radiogenomics/Software/BAKRGibbs.cpp")

######################################################################################
######################################################################################
######################################################################################

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

######################################################################################
######################################################################################
######################################################################################

### Load in the TCGA Gene Expression Data ###
load('~/Dropbox (Personal)/Columbia Radiogenomics/Data/TCGA_GBM_Expression.RData')
G = t(X_TCGA[,which(colnames(X_TCGA)%in%rownames(ECs))])

### Find Overlapping Samples: ECs and Gene Expression ### 
ECs = ECs[which(rownames(ECs)%in%rownames(G)),]
G = G[match(rownames(ECs),rownames(G)),]

######################################################################################
######################################################################################
######################################################################################

### Load in the Morphometric Data For the TCIA MRIs ###
Morph = read.table("~/Dropbox (Personal)/Columbia Radiogenomics/Data/TCIA_Morphometrics.txt",header = TRUE)
rownames(Morph) = paste(Morph$Feature,Morph$Statistics,sep = "_")
colnames(Morph) = gsub("[.]", "-", as.character(colnames(Morph)))
Morph = t(Morph[,-c(1:2)])

### Find Overlapping Samples: ECs and Morphometrics ### 
Morph = Morph[which(rownames(Morph)%in%rownames(ECs)),]
ECs = ECs[which(rownames(ECs)%in%rownames(Morph)),]

### Find Overlapping Samples: Expression and Morphometrics ### 
G = G[which(rownames(G)%in%rownames(Morph)),]

######################################################################################
######################################################################################
######################################################################################

### Load in the Geometric Data For the TCIA MRIs ###
Geo = read.xls("~/Dropbox (Personal)/Columbia Radiogenomics/Data/TCIA_Geometrics.xls")
rownames(Geo) = Geo$Patient.ID; Geo = as.matrix(Geo[,-c(1,6:11)])
Geo = Geo[which(rownames(Geo)%in%rownames(ECs)),]

### Find Overlapping Samples: Geometrics vs. {ECs, Expression, Morphometrics} ### 
ECs = ECs[which(rownames(ECs)%in%rownames(Geo)),]
G = G[which(rownames(G)%in%rownames(Geo)),]
Morph = Morph[which(rownames(Morph)%in%rownames(Geo)),]

### Find Overlapping Samples: ECs vs. {Expression, Morphometrics, Geometrics} ### 
G = G[match(rownames(ECs),rownames(G)),]
Geo = Geo[match(rownames(ECs),rownames(Geo)),]
Morph = Morph[match(rownames(ECs),rownames(Morph)),]

### Check to Make Sure the Sample Sizes are Equal ###
dim(ECs); dim(G); dim(Geo); dim(Morph);

######################################################################################
######################################################################################
######################################################################################

### Read in the Phenotypes ###
Phenos = read.csv("~/Dropbox (Personal)/Columbia Radiogenomics/Data/TCGA_Clinical_Traits.csv")
Y = Phenos; rownames(Y) = as.character(Y$Patient.ID); Y = Y[,-1]

### Find Overlapping Samples: ECs and Clinical Traits (DFS and OS) ### 
Y = Y[which(rownames(Y)%in%rownames(ECs)),17:18]

### Find Overlapping Samples: Clinical Traits vs. {ECs, Expression, Morphometrics, Geometrics} ### 
ECs = ECs[which(rownames(ECs)%in%rownames(Y)),]
G = G[which(rownames(G)%in%rownames(Y)),]
Morph = Morph[which(rownames(Morph)%in%rownames(Y)),]
Geo = Geo[which(rownames(Geo)%in%rownames(Y)),]

### Check the Dimensionalities of the Data ###
dim(ECs); dim(G); dim(Geo); dim(Y); dim(Morph)

######################################################################################
######################################################################################
######################################################################################

### Setup Parameters for Analysis ###

### Set the random seed for reproducibility ###
set.seed(11151990,kind = "L'Ecuyer-CMRG")

#Define the percentage of data to use for splits
nsplit = 0.8; ndatasets = 1e3

#Define the number of datasets to simulate (i.e. times to split)
iter = 2e4
burn = 1e4
thin = 10
sigma=1

### Load in the Cross-Validated Bandwidths ###
load("~/Dropbox (Personal)/Columbia Radiogenomics/Analysis/Cross_Validation_Results/CauchyCV_Results.RData"); 

### Call the Spectrum of Considered Bandwidths ###
theta = seq(from = 0.1, to = 10, by = 0.1)

### Call the available cores accross the computer ###
registerDoParallel(cores = detectCores())

### Set up List ###
Res = list()

######################################################################################
######################################################################################
######################################################################################

### Run the Analysis ###
for(j in 1:ncol(Y)){
  
  ### Choose Phenotype/Trait ###
  y = scale(Y[!is.na(Y[,j]),j])
  
  ### Set up a list to save results ###
  Results = foreach(i = 1:ndatasets)%dopar%{
    
    ### Create the training and test sets ###
    ind = sample(1:length(y), size=nsplit*length(y), replace=FALSE)
    
    ######################################################################################
    ######################################################################################
    ######################################################################################
    
    ### Geometric Data Analysis ###
    X = Geo[!is.na(Y[,j]),]; X = scale(X); X = cbind(rep(1,nrow(X)),X)
    theta_hat = theta[cvs[[j]][2]]
    K = CauchyKernel(t(X),theta_hat); diag(K) = 1
    
    ### Center and Scale the Covariance Matrix ###
    n=nrow(K)
    v=matrix(1, n, 1)
    M=diag(n)-v%*%t(v)/n
    K=M%*%K%*%M
    K=K/mean(diag(K))
    
    ### Compute the Posterior Predictive Mean ###
    Kn = K[ind,ind]
    fhat = K[,ind] %*% solve(Kn + diag(sigma,nrow(Kn)),y[ind])
    
    ### Assess the Model Type with the Mean Absolute Error and Predictive R^2 ###
    MAE.Geo = mean(abs(y[-ind]-fhat[-ind]))
    R2_Geo = cor(y[-ind],fhat[-ind])^2
    
    ######################################################################################
    ######################################################################################
    ######################################################################################
    
    ### Morphometric Data Analysis ###
    X = Morph[!is.na(Y[,j]),]; X = scale(X,scale=FALSE); X = cbind(rep(1,nrow(X)),X)
    theta_hat = theta[cvs[[j]][3]]
    K = CauchyKernel(t(X),theta_hat); diag(K) = 1
    
    ### Center and Scale the Covariance Matrix ###
    n=nrow(K)
    v=matrix(1, n, 1)
    M=diag(n)-v%*%t(v)/n
    K=M%*%K%*%M
    K=K/mean(diag(K))
    
    ### Compute the Posterior Predictive Mean ###
    Kn = K[ind,ind]
    fhat = K[,ind] %*% solve(Kn + diag(sigma,nrow(Kn)),y[ind])
    
    ### Assess the Model Type with the Mean Absolute Error and Predictive R^2 ###
    MAE.Morph = mean(abs(y[-ind]-fhat[-ind]))
    R2_Morph = cor(y[-ind],fhat[-ind])^2
    
    ######################################################################################
    ######################################################################################
    ######################################################################################
    
    ### Gene Expression Analysis ###
    X = G[!is.na(Y[,j]),]; X = log2(X+1); X = cbind(rep(1,nrow(X)),X)
    theta_hat = theta[cvs[[j]][1]]
    K = CauchyKernel(t(X),theta_hat); diag(K) = 1
    
    ### Center and Scale the Covariance Matrix ###
    n=nrow(K)
    v=matrix(1, n, 1)
    M=diag(n)-v%*%t(v)/n
    K=M%*%K%*%M
    K=K/mean(diag(K))
    
    ### Compute the Posterior Predictive Mean ###
    Kn = K[ind,ind]
    fhat = K[,ind] %*% solve(Kn + diag(sigma,nrow(Kn)),y[ind])
    
    ### Assess the Model Type with the Mean Absolute Error and Predictive R^2 ###
    MAE.G = mean(abs(y[-ind]-fhat[-ind]))
    R2_G = cor(y[-ind],fhat[-ind])^2
    
    ######################################################################################
    ######################################################################################
    ######################################################################################
    
    ### Topological Summary Statistics Analysis ###
    X = ECs[!is.na(Y[,j]),]; X = scale(X); X = cbind(rep(1,nrow(X)),X)
    theta_hat = theta[cvs[[j]][4]]
    K = CauchyKernel(t(X),theta_hat); diag(K) = 1
    
    ### Center and Scale the Covariance Matrix ###
    n=nrow(K)
    v=matrix(1, n, 1)
    M=diag(n)-v%*%t(v)/n
    K=M%*%K%*%M
    K=K/mean(diag(K))
    
    ### Compute the Posterior Predictive Mean ###
    Kn = K[ind,ind]
    fhat = K[,ind] %*% solve(Kn + diag(sigma,nrow(Kn)),y[ind])
    
    ### Assess the Model Type with the Mean Absolute Error and Predictive R^2 ###
    MAE.ECs = mean(abs(y[-ind]-fhat[-ind]))
    R2_ECs = cor(y[-ind],fhat[-ind])^2
    
    ######################################################################################
    ######################################################################################
    ######################################################################################
    
    ### Save The Results ###
    c(MAE.G,MAE.Morph,MAE.Geo,MAE.ECs,R2_G,R2_Morph,R2_Geo,R2_ECs)
  }
  
  ### Create a matrix to save results ###
  Final = matrix(unlist(Results),nrow = ndatasets,ncol = 8,byrow = TRUE)
  mod.names = c("Expression","Morph","Geo","ECs")
  colnames(Final) = c(paste("RMSEP",mod.names,sep="_"),paste("R2",mod.names,sep="_"))
  
  ### Save the Results ###
  Res[[j]] = Final; cat("Completed Phenotype",j,"\n")
}

### Check the Results: (1) Overall Survival; (2) Disease Free Surival ###
colMeans(Res[[1]][,1:4]); colMeans(Res[[1]][,5:8])
colMeans(Res[[2]][,1:4]); colMeans(Res[[2]][,5:8])

### Save the Results ###
save(Res,file = "~/Dropbox (Personal)/Columbia Radiogenomics/Analysis/Results/Cauchy_SECT.RData")
