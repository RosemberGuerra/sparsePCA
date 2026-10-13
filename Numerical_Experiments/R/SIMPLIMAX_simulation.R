### Title:    sIMPLIMAX simulation
### Author:   Rosember Guerrag
### Created:  19-09-2019
### Modified: 10-03-2019


rm(list = ls(all.names = TRUE))

library(doParallel)


# Setting the working directory to the source file location
current_working_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(current_working_dir)
getwd()

# Registering the cluster
no_cores <- detectCores() - 1
c1 <- makePSOCKcluster(floor( no_cores*.8))
registerDoParallel(c1)

# Simulation #
# 1 means the model that match the data set   #
# 2 means P and W sparse and equals           #
# 3 means the model and the data don't match  #


DataDir = "S:/DataGeneration"

# Directories - P sparse
mainDir1 = "../../DataGeneration/Psparse/DATA-R"

# Setting directory -P and W sparse 
mainDir2 = "../../DataGeneration/PandWsparse/DATA-R"

# Setting directory - P sparse - Cross data
mainDir3 = "../../DataGeneration/Wsparse/DATA-R"


source("RelativeError.R")
source("Correlation.R")
source("PCAsimplimax.R")


load(paste0(mainDir1,"/Info_simulaiton.RData"))

Ndatasets = Infor_simulation$n_data_sets
MatrixInfo = Infor_simulation$design_matrix_replication

Simulation = foreach(i=1:Ndatasets,
                     .packages = c("combinat","psych","GPArotation"),
                     .combine=rbind)%dopar%{
                       
                       # Loading data
                       setwd(file.path(DataDir))
                       load(paste0("Psparse/DATA-R/Psparse",i,".RData")) # P sparse
                       Out1 = list(X = out$X, P =out$P, W = out$W, Z= out$Z, 
                                   k = out$k,Propsparse= out$Propsparse )
                       
                       load(paste0("PandWsparse/DATA-R/PWsparse",i,".RData")) # P and W sparse
                       Out2 = list(X = out$X, P =out$P, W = out$W, Z= out$Z, 
                                   k = out$k,Propsparse= out$Propsparse )
                       
                       
                       load(paste0("Wsparse/DATA-R/Wsparse",i,".RData")) # W sparse
                       Out3 = list(X = out$X, P =out$P, W = out$W, Z= out$Z,
                                   k = out$k,Propsparse= out$Propsparse )
                       
                       
                       PZ1 = PCAsimplimax(X = Out1$X,k = Out1$k , ps = Out1$Propsparse )   # P sparse
                       PZ2 = PCAsimplimax(X = Out2$X,k = Out2$k , ps = Out2$Propsparse )   # P and W sparse
                       PZ3 = PCAsimplimax(X = Out3$X,k = Out3$k , ps = Out3$Propsparse )   # W sparse
                       
                       P1 = PZ1$P
                       P2 = PZ2$P
                       P3 = PZ3$P
                       
                       Z1 = PZ1$Z
                       Z2 = PZ2$Z
                       Z3 = PZ3$Z
                       
                       out_error1 =  RelativeError(Out1$P, P1) # P sparse
                       out_error2 =  RelativeError(Out2$P, P2) # P and W sparse
                       out_error3 =  RelativeError(Out3$P, P3) # W sparse- Cross Data
                       
                       P1 = P1[,out_error1$permutation]%*%diag(out_error1$s)
                       P2 = P2[,out_error2$permutation]%*%diag(out_error2$s)
                       P3 = P3[,out_error3$permutation]%*%diag(out_error3$s)
                       
                       Z1 = Z1[,out_error1$permutation]%*%diag(out_error1$s)
                       Z2 = Z2[,out_error2$permutation]%*%diag(out_error2$s)
                       Z3 = Z3[,out_error3$permutation]%*%diag(out_error3$s)
                       
                       ErrorLW1 = out_error1$Error
                       ErrorLW2 = out_error2$Error
                       Correlation3 = correlation(Out3$W, P3)  
                       
                       
                       ErrorScores1 = (norm((Out1$Z-Z1),type = "F")/norm(Out1$Z,type = "F"))^2
                       ErrorScores2 = (norm((Out2$Z-Z2),type = "F")/norm(Out2$Z,type = "F"))^2
                       CorrelationZ3 = correlation(Out3$Z, Z3) 
                       
                       ZeroInd1 = which(as.vector(Out1$P)==0) # P sparse 
                       ZeroInd2 = which(as.vector(Out2$P)==0) # P sparse 
                       ZeroInd3 = which(as.vector(Out3$W)==0) # W sparse 
                       
                       W1_vec = as.vector(P1)
                       W2_vec = as.vector(P2)
                       W3_vec = as.vector(P3)
                       
                       if (Out1$Propsparse == 0){ # P sparse
                         MR01 = 0
                         MR03 = 0
                       }else{
                         wz1 = W1_vec[ZeroInd1]
                         wz3 = W3_vec[ZeroInd3]
                         
                         MR01 = 1-sum(wz1 == 0)/length(ZeroInd1)
                         MR03 = 1-sum(wz3 == 0)/length(ZeroInd3)
                       }
                       if (Out2$Propsparse == 0){ # P sparse
                         MR02 = 0
                       }else{
                         wz2 = W2_vec[ZeroInd2]
                         MR02 = 1-sum(wz2 == 0)/length(ZeroInd2)
                       }
                       
                       
                       
                       
                       # Sparse Variance #
                       Xs1 = Z1%*%t(P1)
                       Xs2 = Z2%*%t(P2)
                       Xs3 = Z3%*%t(P3)
                       
                       Svar1 = 1 - (norm((Out1$X-Xs1),type = "F")/norm(Out1$X,type = "F"))^2
                       Svar2 = 1 - (norm((Out2$X-Xs2),type = "F")/norm(Out2$X,type = "F"))^2
                       Svar3 = 1 - (norm((Out3$X-Xs3),type = "F")/norm(Out3$X,type = "F"))^2
                       
                       
                       # Saving performance #
                       return(list(ErrorLW1 = ErrorLW1, ErrorScores1 = ErrorScores1, MR01 = MR01,Svar1 = Svar1,
                                   ErrorLW2 = ErrorLW2, ErrorScores2 = ErrorScores2, MR02 = MR02,Svar2 = Svar2,
                                   Correlation3= Correlation3, CorrelationZ3= CorrelationZ3, MR03 = MR03,Svar3 = Svar3))
                       
                     }

# stop Cluster
stopCluster(c1)

## Unlisting the results ##

ErrorLW1 = as.numeric(unlist(Simulation[,"ErrorLW1"]))
ErrorLW2 = as.numeric(unlist(Simulation[,"ErrorLW2"]))
Correlation3 = as.numeric(unlist(Simulation[,"Correlation3"]))

ErrorScores1 = as.numeric(unlist(Simulation[,"ErrorScores1"]))
ErrorScores2 = as.numeric(unlist(Simulation[,"ErrorScores2"]))
CorrelationZ3 = as.numeric(unlist(Simulation[,"CorrelationZ3"]))


MR01 = as.numeric(unlist(Simulation[,"MR01"]))
MR02 = as.numeric(unlist(Simulation[,"MR02"]))
MR03 = as.numeric(unlist(Simulation[,"MR03"]))


Svar1 = as.numeric(unlist(Simulation[,"Svar1"]))
Svar2 = as.numeric(unlist(Simulation[,"Svar2"]))
Svar3 = as.numeric(unlist(Simulation[,"Svar3"]))

## Results in a table ##
Methodology = rep('Simplimax',Ndatasets)

load("../../DataGeneration/PandWsparse/DATA-R/Info_simulaiton.RData")
p_sparse2 = Infor_simulation$design_matrix_replication[,"p_sparse"]

DataSimplimax = cbind(Methodology,p_sparse2, MatrixInfo,
                 ErrorLW1,ErrorScores1,MR01, Svar1,
                 ErrorLW2,ErrorScores2,MR02, Svar2,
                 Correlation3,CorrelationZ3,MR03, Svar3)
write.table(DataSimplimax, file = "JointDataSimplimax.txt", sep = ",",row.names = FALSE)
save.image(file = "SimplimaxResults.RData")
