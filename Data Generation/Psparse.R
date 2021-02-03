# P sparse using foreach
# Rosember Guerra
# 15-09-2019

# install.packages("doParallel")    # Install doParallel package
# install.packages("MASS")          # Install MASS package
# install.packages("mvtnorm")       # Install mvtnorm package

rm(list = ls(all.names = TRUE))
current_working_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(current_working_dir)
getwd()

dir.create("DATA-R") # Directory to save the data

# setting the number of cores 
library(doParallel)

no_cores <- detectCores() 
c1 <- makePSOCKcluster(floor( no_cores*.8))
registerDoParallel(c1)

# sparse PCA data simulation #
# sparse PCA data simulation #

set.seed(2019)

# sparse PCA data simulation #
VAFx = c(.80,.95,1)           # Proportion of explained variance
p_sparse = c(0,.5,.8)      # Proportion of sparsity
n_components = c(2,3)         # Number of components
s_size = c(100,500)           # Sample size
n_variables = c(10,100,1000)     # Number of variables
n_replications = c(100)        # Number of repetitions



design_matrix <- expand.grid( n_variables=n_variables,s_size=s_size,p_sparse=p_sparse,
                              n_components=n_components,VAFx =VAFx)

design_matrix_replication <- design_matrix[rep(1:nrow(design_matrix), times = n_replications), ]

Infor_simulation = list(n_data_sets = nrow(design_matrix_replication), n_replications  =n_replications,
                        design_matrix_replication = design_matrix_replication)
save(Infor_simulation, file = "DATA-R/Info_simulaiton.RData")


# start simulating the data 

results_sim1_data1 <- foreach(i=1:nrow(design_matrix_replication),
                              .options.RNG = 2008,
                              .packages = c("MASS","mvtnorm"),
                              .combine=rbind)%dopar%{
                                
                                # set the specific values of the parameters
                                R = design_matrix_replication$n_components[i]
                                I =  design_matrix_replication$s_size[i]
                                J = design_matrix_replication$n_variables[i]
                                vafx = design_matrix_replication$VAFx[i]
                                propsparse = design_matrix_replication$p_sparse[i]
                                
                                n_zeros = floor(J*propsparse)
                                n_non_zero = J-n_zeros
                                
                                # Generating data
                                X = mvrnorm(n =I,mu=rep(0,J),Sigma = diag(1,J))
                                
                                # SVD on Xinit
                                svd1  =  svd(X)
                                
                                # P matrix is the right singular vectors
                                P = svd1$v[,1:R]
                                if (propsparse !=0){
                                  P = apply(P, 2, function(x){
                                    t = quantile(abs(x),propsparse)
                                    x[abs(x) < t] = 0
                                    return(x)
                                  })
                                }
                                
                                #
                                
                                TruecScores = svd1$u[,1:R]
                                
                                #
                                P = P%*%diag(svd1$d[1:R])
                                Xtrue = TruecScores%*%t(P)
                                W = ginv(Xtrue)%*%TruecScores
                                 # Adding noice
                                SSqXtrue =  sum(Xtrue^2)                        # sum squares of the data set
                                EX = mvrnorm(I,mu= rep(0,J),Sigma = diag(1,J))  # EX = Error of X
                                SSqEX = sum(EX^2)                               # Sum squares fo the EX
                                fx = sqrt(SSqXtrue*(1-vafx)/(vafx * SSqEX))
                                Xnew = Xtrue + fx*EX                            # Data with noise
                                
                                # Saving data
                                out = list(X = Xnew, W =W ,P = P, Z = TruecScores, k = R, 
                                           Propsparse = propsparse, nzeros =  n_zeros, teller = i)
                                save(out, file=paste0("DATA-R/Psparse",i, ".RData"))
                                
                              }

# stop Cluster
stopCluster(c1)
