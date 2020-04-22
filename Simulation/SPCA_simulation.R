### Title:    SPCA simulation 
### Author:   Rosember Guerrag
### Created:  19-09-2019
### Modified: 10-03-2019


# rm(list = ls(all.names = TRUE))

library(doParallel)

# Setting the working directory to the source file location
current_working_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(current_working_dir)
getwd()

# Registering the cluster
no_cores <- detectCores()
c1 <- makePSOCKcluster(floor( no_cores*.8))
registerDoParallel(c1)

# The function trace is used for not standarized the data set as follow #
# trace(elasticnet::spca, edit = TRUE) #

# Simulation #
# 1 means the model that match the data set   #
# 2 means P and W sparse and equals           #
# 3 means the model and the data don't match  # 

# Setting directory - W sparse 
DataDir = "S:/DataGeneration"

mainDir1 = "../../DataGeneration/Wsparse/DATA-R"
#ResultDir1 = "../Simulation/SPCA/Results/SPCA1/"

# Setting directory -P and W sparse 
mainDir2 = "../../DataGeneration/PandWsparse/DATA-R"
#ResultDir2 = "../Simulation/SPCA/Results/SPCA2/"

# Setting directory - P sparse - Cross data
mainDir3 = "../../DataGeneration/Psparse/DATA-R"
#ResultDir3 = "../Simulation/SPCA/Results/SPCA3/"


load(paste0(mainDir1,"/Info_simulaiton.RData"))
Info_matrix = Infor_simulation$design_matrix_replication

source("RelativeError.R")
source("Correlation.R")

Ndatasets = Infor_simulation$n_data_sets
MatrixInfo = Infor_simulation$design_matrix_replication

Simulation = foreach(i=1:Ndatasets,.packages = c("elasticnet","combinat","MASS"),
        .combine=rbind)%dopar%{
          
          # Loading data
          setwd(file.path(DataDir))
          load(paste0("Wsparse/DATA-R/Wsparse",i,".RData")) # W sparse
          Out1 = list(X = out$X, P =out$P, W = out$W, Z= out$Z, nzeros = out$nzeros,
                      k = out$k,Propsparse= out$Propsparse )
          
          load(paste0("PandWsparse/DATA-R/PWsparse",i,".RData")) # P and W sparse
          Out2 = list(X = out$X2, P =out$P, W = out$W, Z= out$Z2, nzeros = out$nzeros,
                      k = out$k,Propsparse= out$Propsparse )
          
          load(paste0("Psparse/DATA-R/Psparse",i,".RData")) # P sparse
          Out3 = list(X = out$X, P =out$P, W = out$W, Z= out$Z, nzeros = out$nzeros,
                      k = out$k,Propsparse= out$Propsparse )
          
          # Estimation of SPCA #
          n_non_zeros1 = dim(Out1$X)[2]-Out1$nzeros
          n_non_zeros2 = dim(Out2$X)[2]-Out2$nzeros
          n_non_zeros3 = dim(Out3$X)[2]-Out3$nzeros
          
          Wmatrix1 =  spca(x=Out1$X,K = Out1$k, sparse = "varnum", para =rep(n_non_zeros1,Out1$k) ,type = "predictor")
          Wmatrix2 =  spca(x=Out2$X,K = Out2$k, sparse = "varnum", para =rep(n_non_zeros2,Out2$k) ,type = "predictor")
          Wmatrix3 =  spca(x=Out3$X,K = Out3$k, sparse = "varnum", para =rep(n_non_zeros3,Out3$k) ,type = "predictor")
          
          W1 = Wmatrix1$loadings
          W2 = Wmatrix2$loadings
          W3 = Wmatrix3$loadings
          
          Z1 = Out1$X%*%W1
          Z2 = Out2$X%*%W2
          Z3 = Out3$X%*%W3
          
          out_error1 =  RelativeError(Out1$W, W1) # W sparse
          out_error2 =  RelativeError(Out2$W, W2) # P and W sparse
          out_error3 =  RelativeError(Out3$W, W3) # P sparse- Cross Data
          
          W1 = W1[,out_error1$permutation]%*%diag(out_error1$s)
          W2 = W2[,out_error2$permutation]%*%diag(out_error2$s)
          W3 = W3[,out_error3$permutation]%*%diag(out_error3$s)
          
          Z1 = Z1[,out_error1$permutation]%*%diag(out_error1$s)
          Z2 = Z2[,out_error2$permutation]%*%diag(out_error2$s)
          Z3 = Z3[,out_error3$permutation]%*%diag(out_error3$s)
          
          ErrorLW1 = out_error1$Error
          ErrorLW2 = out_error2$Error
          #ErrorLW3 = out_error3$Error
          Correlation3 = correlation(Out3$P, W3) 
          
          ErrorScores1 = (norm((Out1$Z-Z1),type = "F")/norm(Out1$Z,type = "F"))^2
          ErrorScores2 = (norm((Out2$Z-Z2),type = "F")/norm(Out2$Z,type = "F"))^2
          #ErrorScores3 = (norm((Out3$Z-Z3),type = "F")/norm(Out3$Z,type = "F"))^2
          CorrelationZ3 = correlation(Out3$Z, Z3) 
          
          ZeroInd1 = which(as.vector(Out1$W)==0) # W sparse 
          ZeroInd2 = which(as.vector(Out2$W)==0) # W sparse 
          ZeroInd3 = which(as.vector(Out3$P)==0) # P sparse 
          
          W1_vec = as.vector(W1)
          W2_vec = as.vector(W2)
          W3_vec = as.vector(W3)
          
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
          P1 = t(Out1$X)%*%ginv(t(Z1))
          P2 = t(Out2$X)%*%ginv(t(Z2))
          P3 = t(Out3$X)%*%ginv(t(Z3))
          
          Xs1 = Z1%*%t(P1)
          Xs2 = Z2%*%t(P2)
          Xs3 = Z3%*%t(P3)
          
          Svar1 = 1 - (norm((Out1$X-Xs1),type = "F")/norm(Out1$X,type = "F"))^2
          Svar2 = 1 - (norm((Out2$X-Xs2),type = "F")/norm(Out2$X,type = "F"))^2
          Svar3 = 1 - (norm((Out3$X-Xs3),type = "F")/norm(Out3$X,type = "F"))^2
          
          # saving results
          # results1 = list(P = P1, W = W1, Z = Z1, PEV =PEV1)
          # results2 = list(P = P2, W = W2, Z = Z2, PEV =PEV2)
          # results3 = list(P = P3, W = W3, Z = Z3, PEV =PEV3)
          # save(results1, file = paste0(ResultDir1,"Results",i,".RData"))
          # save(results2, file = paste0(ResultDir2,"Results",i,".RData"))
          # save(results3, file = paste0(ResultDir3,"Results",i,".RData"))
          
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
#ErrorLW3 = as.numeric(unlist(Simulation[,"ErrorLW3"]))
Correlation3 = as.numeric(unlist(Simulation[,"Correlation3"]))

ErrorScores1 = as.numeric(unlist(Simulation[,"ErrorScores1"]))
ErrorScores2 = as.numeric(unlist(Simulation[,"ErrorScores2"]))
#ErrorScores3 = as.numeric(unlist(Simulation[,"ErrorScores3"]))
CorrelationZ3 = as.numeric(unlist(Simulation[,"CorrelationZ3"]))

MR01 = as.numeric(unlist(Simulation[,"MR01"]))
MR02 = as.numeric(unlist(Simulation[,"MR02"]))
MR03 = as.numeric(unlist(Simulation[,"MR03"]))


Svar1 = as.numeric(unlist(Simulation[,"Svar1"]))
Svar2 = as.numeric(unlist(Simulation[,"Svar2"]))
Svar3 = as.numeric(unlist(Simulation[,"Svar3"]))

## Results in a table ##
Methodology = rep('SPCA',Ndatasets)

load("../../DataGeneration/PandWsparse/DATA-R/Info_simulaiton.RData")
p_sparse2 = Infor_simulation$design_matrix_replication[,"p_sparse"]

DataSPCA = cbind(Methodology,p_sparse2, MatrixInfo,
                 ErrorLW1,ErrorScores1,MR01, Svar1,
                 ErrorLW2,ErrorScores2,MR02, Svar2,
                 Correlation3,CorrelationZ3,MR03, Svar3)

write.table(DataSPCA, file = "JointDataSPCA2.txt", sep = ",",row.names = FALSE)
save.image(file = "SPCAresults.RData")
