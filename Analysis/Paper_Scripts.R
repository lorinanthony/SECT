### The purpose of this file is to reproduce the predictive GBM survival analysis and main results presented in Table 1 from the SECT manuscript. Before running this script, please make sure that your current working directory contains all elements from the downloaded SECT repo ###

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

### NOTE: Please run each of these lines in subsequent order. The last lines should print a table. Here, we predict disease free survival (DFS) and overall survival (OS) using Gaussian process (GP) regression models defined by the linear, Gaussian, and Cauchy covariance functions, respectively. For each model fit, we consider the predictive utility of four different genomic data types: gene expression, tumor morphometry, tumor geometry, and the proposed smooth Euler characteristic transform (SECT). Assessment is carried out by using the predictive squared correlation coefficient (R^2), where larger numbers indicate better performance. We also use Optimal% to denote the percentage of the time that a model exhibits the greatest R^2. All values in bold represent the best method in these two assessment categories. These values are based on 1000 random 80-20 splits for each clinical outcome. Standard errors for each model are given the parentheses. Lastly, we give estimates for the bandwidth or length-scale parameter θ used to compute each kernel function. Note that θ was found by using 10-fold cross-validation over the grid [0.1, 10] with step sizes equal to 0.1. ###

### GP Regression w/ Linear Covariance Function ###
source("./Analysis/Linear_Function_Analysis.R")

### GP Regression w/ Gaussian Covariance Function ###
source("./Analysis/Gauss_Function_Analysis.R")

### GP Regression w/ Cauchy Covariance Function ###
source("./Analysis/Cauchy_Function_Analysis.R")

######################################################################################
######################################################################################
######################################################################################

### Present Prediction Results ###
source("./Analysis/Prediction_Results/Results.R")
