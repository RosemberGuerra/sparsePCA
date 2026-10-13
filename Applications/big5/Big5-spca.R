### Title:    Application Big5 data set.
### Model:    SPCA 
### Author:   Rosember Guerra-Urzola
### Created:  16-11-2020
### Modified: 

library(qgraph)
library(psych)
library(latticeExtra)
library(plotly)
library(elasticnet)
library(MASS)
current_working_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(current_working_dir)
getwd()

rm(list = ls(all.names = TRUE))

#functions#
# source("PCAvarimax.R")
# source("PCAsimplimax.R")

set.seed(1234)
# loading data set #
data("big5")
DIM = dim(big5)
X = scale(big5, center = TRUE, scale = TRUE ) 


Xsvd = svd(X)
# plot - pev vs components $
pev = (Xsvd$d^2)/sum(Xsvd$d^2)
plot(cumsum(pev),type = "l",ylab = "% of explained variance", xlab = "Number of components",
     main = "Big 5")

K = 5 # number of components
Xpca = Xsvd$u[,1:K]%*%diag(Xsvd$d[1:K])%*%t(Xsvd$v[,1:K])
Rv0 = 1 - (norm((Xpca-X),type = "F")/norm(X,type = "F"))^2

Wmatrix1=spca(x=X,K=K,sparse = "varnum",para =rep(64,K),
              type = "predictor")
Z=X%*%Wmatrix1$loadings
P=t(X)%*%ginv(t(Z))
Xs = Z%*%t(P)
pevspca= 1-norm((Xs-X),type = "F")/norm(X,type = "F")
W=Wmatrix1$loadings

# table #
O=grep('O',colnames(X))
C=grep('C',colnames(X)) 
E=grep('E',colnames(X))
A=grep('A',colnames(X))
N=grep('N',colnames(X))

nonzeroOCEAN=function(P){
  return( c(sum(P[O]!=0),sum(P[C]!=0),sum(P[E]!=0),sum(P[A]!=0),sum(P[N]!=0)))
}

spca_ocean= apply(W, 2, nonzeroOCEAN)
colSums(spca_ocean)
View(spca_ocean)
