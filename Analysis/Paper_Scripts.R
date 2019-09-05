### The purpose of this file is to reproduce the predictive GBM survival analysis and main results presented in Table 1 from the SECT manuscript. Before running this script, please make sure that your current working directory contains all elements from the downloaded SECT repo ###

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

######################################################################################
######################################################################################
######################################################################################

### NOTE: Please run each of these lines in subsequent order. The last lines should print a table. Here, we predict disease free survival (DFS) and overall survival (OS) using Gaussian process (GP) regression models defined by the linear, Gaussian, and Cauchy covariance functions, respectively. For each model fit, we consider the predictive utility of four different genomic data types: gene expression, tumor morphometry, tumor geometry, and the proposed smooth Euler characteristic transform (SECT). Assessment is carried out in the first column by using the predictive squared correlation coefficient (R^2), where larger numbers indicate better performance. These values are based on 1000 random 80-20 splits for each clinical outcome. Standard errors for each model are given the second column. In the third column, we use Optimal% to denote the percentage of the time that a model exhibits the greatest R^2. Lastly, we give estimates for the bandwidth (or length-scale) parameter θ used to compute each kernel function in the fourth column. Note that θ was found by using 10-fold cross-validation over the grid [0.1, 10] with step sizes equal to 0.1. ###

### Update: We have added two new lines to track how long it takes to run the analyses. It should take about 3.5 minutes on a MacBook Pro (Processor: 3.1-gigahertz (GHz) Intel Core i5, Memory: 8GB 2133-megahertz (MHz) LPDDR3).

ptm <- proc.time() #Start clock for timing

### GP Regression w/ Linear Covariance Function ###
source("./Analysis/Linear_Function_Analysis.R")

### GP Regression w/ Gaussian Covariance Function ###
source("./Analysis/Gauss_Function_Analysis.R")

### GP Regression w/ Cauchy Covariance Function ###
source("./Analysis/Cauchy_Function_Analysis.R")

proc.time() - ptm #Stop clock and check

######################################################################################
######################################################################################
######################################################################################

### Present Prediction Results ###
source("./Analysis/Prediction_Results/Results.R")
