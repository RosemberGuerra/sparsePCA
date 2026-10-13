### Title:    Application Big5 data set.
### Model:    Varimax and Simplimax
### Author:   Rosember Guerra-Urzola
### Created:  22-01-2020
### Modified: 13-11-2020

library(qgraph)
library(psych)
library(latticeExtra)
library(plotly)


current_working_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(current_working_dir)
getwd()

rm(list = ls(all.names = TRUE))

#functions#
source("PCAvarimax.R")
source("PCAsimplimax.R")

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

VarMx = PCAvarimax(X,K,.734)
Xvmx = VarMx$Z%*%t(VarMx$P)
Rva_varimax = 1- (norm((Xvmx-X),type = "F")/norm(X,type = "F"))^2

SimMx = PCAsimplimax(X,K,.734)
Xspx = SimMx$Z%*%t(SimMx$P)
Rva_simplimax = 1- (norm((Xspx-X),type = "F")/norm(X,type = "F"))^2

save(VarMx,Xvmx,file = "VarimaxApplication.RData")
save(SimMx,Xspx,file = "SimplimaxApplication.RData")

# table #
O=grep('O',colnames(X))
C=grep('C',colnames(X)) 
E=grep('E',colnames(X))
A=grep('A',colnames(X))
N=grep('N',colnames(X))

nonzeroOCEAN=function(P){
  return( c(sum(P[O]!=0),sum(P[C]!=0),sum(P[E]!=0),sum(P[A]!=0),sum(P[N]!=0)))
}

V_ocean= apply(VarMx$P, 2, nonzeroOCEAN)
colSums(V_ocean)
View(V_ocean)
S_ocean= apply(SimMx$P, 2, nonzeroOCEAN)
colSums(S_ocean)
View(S_ocean)

## correlation between varimax and simplimax ##
