# Joint data #
# Rosember Guerra 
# update: 02-01-2020

library(data.table)
library(ggplot2)
library(ggpubr)

current_working_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(current_working_dir)
getwd()

rm(list = ls(all.names = TRUE))

# loading the data sets #
# These data are the final files that comen from the simulation in each methodology

Simplimax = read.csv("JointDataSimplimax.txt", sep = ",", header = TRUE)
Varimax = read.csv("JointDataVarimax.txt", sep = ",", header = TRUE)
sPCArSVD = read.csv("JointDatasPCArSVD.txt", sep = ",", header = TRUE)
#
SPCA = read.csv("JointDataSPCA2.txt", sep = ",", header = TRUE)
pathSPCA = read.csv("JointDataPathSPCA.txt", sep = ",", header = TRUE)
GPower = read.csv("JointDataGPower.txt", sep = ",", header = TRUE)
#

# Join the data #

JoinDataSet = rbind(Simplimax,Varimax,sPCArSVD,SPCA,pathSPCA,GPower)

JoinDataSet = setDT(JoinDataSet)
class(JoinDataSet)

# factor variables #
JoinDataSet$Methodology = as.factor(JoinDataSet$Methodology)
levels(JoinDataSet$Methodology)

JoinDataSet$p_sparse2 = as.factor(JoinDataSet$p_sparse2)
levels(JoinDataSet$p_sparse2)
levels(JoinDataSet$p_sparse2) = c("Sparsity 70%","Sparsity 80%","Sparsity 90%")
levels(JoinDataSet$p_sparse2)

JoinDataSet$n_variables = as.factor(JoinDataSet$n_variables)
levels(JoinDataSet$n_variables)
levels(JoinDataSet$n_variables) = c("J = 10","J = 100", "J = 1000")
levels(JoinDataSet$n_variables)

JoinDataSet$s_size = as.factor(JoinDataSet$s_size)
levels(JoinDataSet$s_size)
levels(JoinDataSet$s_size) = c("I = 100", "I = 500")
levels(JoinDataSet$s_size)

JoinDataSet$p_sparse = as.factor(JoinDataSet$p_sparse)
levels(JoinDataSet$p_sparse)
levels(JoinDataSet$p_sparse) = c("Sparsity 0%","Sparsity 50%","Sparsity 80%")
levels(JoinDataSet$p_sparse)

JoinDataSet$n_components = as.factor(JoinDataSet$n_components)
levels(JoinDataSet$n_components)
levels(JoinDataSet$n_components) = c("K = 2", "K = 3")
levels(JoinDataSet$n_components)

JoinDataSet$VAFx = as.factor(JoinDataSet$VAFx)
levels(JoinDataSet$VAFx)
levels(JoinDataSet$VAFx) = c("vafx = 80%","vafx = 95%","vafx = 100%" )
levels(JoinDataSet$VAFx)


# Second data set # 

LevesFactors = JoinDataSet[rep(1:nrow(JoinDataSet),times =4),1:7]
Measurements = as.factor(as.vector(cbind(rep("SRE-LW",nrow(JoinDataSet)),
                    rep("SRE-S",nrow(JoinDataSet)),
                    rep("PEV",nrow(JoinDataSet)),
                    rep("MR",nrow(JoinDataSet)))))
Measurements = factor(Measurements, levels(Measurements)[c(3,4,2,1)])
levels(Measurements)

Measurements2 = as.factor(as.vector(cbind(rep("CosSim-LW",nrow(JoinDataSet)),
                                          rep("CosSim-S",nrow(JoinDataSet)),
                                          rep("PEV",nrow(JoinDataSet)),
                                          rep("MR",nrow(JoinDataSet)))))
Measurements2 = factor(Measurements2, levels(Measurements2)[c(1,2,4,3)])
levels(Measurements2)


ExpI= as.numeric(as.vector(cbind(JoinDataSet$ErrorLW1,
                                 JoinDataSet$ErrorScores1,
                                 JoinDataSet$Svar1,
                                 JoinDataSet$MR01)))
ExpII = as.numeric(as.vector(cbind(JoinDataSet$ErrorLW2,
                                 JoinDataSet$ErrorScores2,
                                 JoinDataSet$Svar2,
                                 JoinDataSet$MR02)))
ExpIII = as.numeric(as.vector(cbind(JoinDataSet$Correlation3,
                                 JoinDataSet$CorrelationZ3,
                                 JoinDataSet$Svar3,
                                 JoinDataSet$MR03)))

JoinDataSet2 = cbind(LevesFactors,Measurements,Measurements2,ExpI,ExpII,ExpIII)


save(JoinDataSet2, file = "JoinDataSet2.RData")
