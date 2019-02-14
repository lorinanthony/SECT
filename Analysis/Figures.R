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
setwd("~/Dropbox (Personal)/Columbia Radiogenomics/Analysis/Results/")

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
setwd("~/Dropbox (Personal)/Columbia Radiogenomics/Analysis/Figures")

### Look at the Squared Correlation Coefficient ###
dfs.cor = c(colMeans(A[[2]][,5:8]),colMeans(B[[2]][,5:8]),colMeans(C[[2]][,5:8]))
os.cor = c(colMeans(A[[1]][,5:8]),colMeans(B[[1]][,5:8]),colMeans(C[[1]][,5:8]))

dfs.cor.sd = c(apply(A[[2]][,5:8],2,sd),apply(B[[2]][,5:8],2,sd),apply(C[[2]][,5:8],2,sd))/10
os.cor.sd = c(apply(A[[1]][,5:8],2,sd),apply(B[[1]][,5:8],2,sd),apply(C[[1]][,5:8],2,sd))/10

### Disease Free Survival (DFS) ###
pdf("DFS_COR.pdf")
par(mar=c(5,5,4,2))
barx = barplot(dfs.cor, beside=TRUE, space=c(rep(0,4),1,rep(0,3),1,rep(0,3)),col=rep(c("purple3","blue","forest green","dark red"),4),ylim=c(0,0.5), axis.lty=1, xlab=expression("Covariance Function"~(sigma)), ylab=bquote("Squared Correlation Coefficient ("~R^2~")"),border = FALSE,xaxt = "n")
lines(x = rep(4.5,101),y = seq(0,0.5,0.005),col = "grey60",lty = 2,lwd = 1.5)
lines(x = rep(9.5,101),y = seq(0,0.5,0.005),col = "grey60",lty = 2,lwd = 1.5)
error.bar(barx,dfs.cor,dfs.cor.sd)
axis(1,at = c(2,7,12),labels = c("Linear","Gaussian","Cauchy"))
legend("topleft",legend = c("Gene Expression","Morphometrics","Geometrics","SECT"),col = c("purple3","blue","forest green","dark red"),pch=15,bty = "n",cex=1.2,pt.cex = 1.5)
dev.off()

### Overall Survival ###
pdf("OS_COR.pdf")
par(mar=c(5,5,4,2))
barx = barplot(os.cor, beside=TRUE, space=c(rep(0,4),1,rep(0,3),1,rep(0,3)),col=rep(c("purple3","blue","forest green","dark red"),4),ylim=c(0,0.5), axis.lty=1, xlab=expression("Covariance Function"~(sigma)), ylab=bquote("Squared Correlation Coefficient ("~R^2~")"),border = FALSE,xaxt = "n")
lines(x = rep(4.5,101),y = seq(0,0.5,0.005),col = "grey60",lty = 2,lwd = 1.5)
lines(x = rep(9.5,101),y = seq(0,0.5,0.005),col = "grey60",lty = 2,lwd = 1.5)
error.bar(barx,os.cor,os.cor.sd)
axis(1,at = c(2,7,12),labels = c("Linear","Gaussian","Cauchy"))
legend("topleft",legend = c("Gene Expression","Morphometrics","Geometrics","SECT"),col = c("purple3","blue","forest green","dark red"),pch=15,bty = "n",cex=1.2,pt.cex = 1.5)
dev.off()

######################################################################################
######################################################################################
######################################################################################

### Run T-Test to Assess Significance Between Runs ###

### Disease Free Survival ###
t.test(A[[2]][,8],A[[2]][,5],alternative="greater")$p.value
t.test(A[[2]][,8],A[[2]][,6],alternative="greater")$p.value
t.test(A[[2]][,8],A[[2]][,7],alternative="greater")$p.value

t.test(B[[2]][,8],B[[2]][,5],alternative="greater")$p.value
t.test(B[[2]][,8],B[[2]][,6],alternative="greater")$p.value
t.test(B[[2]][,8],B[[2]][,7],alternative="greater")$p.value

t.test(C[[2]][,8],C[[2]][,5],alternative="greater")$p.value
t.test(C[[2]][,8],C[[2]][,6],alternative="greater")$p.value
t.test(C[[2]][,8],C[[2]][,7],alternative="greater")$p.value

### Overall Survival ###
t.test(A[[1]][,8],A[[1]][,5],alternative="greater")$p.value
t.test(A[[1]][,8],A[[1]][,6],alternative="greater")$p.value
t.test(A[[1]][,8],A[[1]][,7],alternative="greater")$p.value

t.test(B[[1]][,8],B[[1]][,5],alternative="greater")$p.value
t.test(B[[1]][,8],B[[1]][,6],alternative="greater")$p.value
t.test(B[[1]][,8],B[[1]][,7],alternative="greater")$p.value

t.test(C[[1]][,8],C[[1]][,5],alternative="greater")$p.value
t.test(C[[1]][,8],C[[1]][,6],alternative="greater")$p.value
t.test(C[[1]][,8],C[[1]][,7],alternative="greater")$p.value
