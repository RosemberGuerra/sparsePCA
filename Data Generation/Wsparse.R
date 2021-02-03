# W sparse using foreach
# Rosember Guerra
# 15-09-2019

# install.packages("doParallel")    # Install doParallel package
# install.packages("MASS")          # Install MASS package
# install.packages("mvtnorm")       # Install mvtnorm package

rm(list = ls(all.names = TRUE))

current_working_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(current_working_dir)
getwd()

# setting the number of cores 
library(doParallel)
no_cores <- detectCores() 
c1 <- makePSOCKcluster(floor( no_cores*.8))
registerDoParallel(c1)

dir.create("DATA-R") # Directory to save the data

# sparse PCA data simulation #

set.seed(2019)

# sparse PCA data simulation #
VAFx = c(0.8,.95,1)           # Proportion of explained variance
p_sparse = c(0,.5,.8)         # Proportion of sparsity
n_components = c(2,3)         # Number of components
s_size = c(100,500)           # Sample size
n_variables = c(10,100,1000)  # Number of variables
n_replications = c(100)         # Number of repetitions



design_matrix <- expand.grid( n_variables=n_variables,s_size=s_size,p_sparse=p_sparse,
                              n_components=n_components,VAFx =VAFx)

design_matrix_replication <- design_matrix[rep(1:nrow(design_matrix), times = n_replications), ]

Infor_simulation = list(n_data_sets = nrow(design_matrix_replication), n_replications  =n_replications,
                        design_matrix_replication = design_matrix_replication)
save(Infor_simulation, file = "DATA-R/Info_simulaiton.RData")


results_sim1_data1 <- foreach(i=1:nrow(design_matrix_replication),
                              .options.RNG = 2018,
                              .packages = c("MASS","mvtnorm"),
                              .combine=rbind)%dopar%{
                                
                                # set the specific values of the parameters
                                R = design_matrix_replication$n_components[i]
                                I =  design_matrix_replication$s_size[i]
                                J = design_matrix_replication$n_variables[i]
                                vafx = design_matrix_replication$VAFx[i]
                                propsparse = design_matrix_replication$p_sparse[i]
                                
                                n_zeros = floor(J*propsparse)
                                n_non_zero = J -n_zeros
                                
                                # 1. Random sample form a multivariate normal distribution
                                X = mvrnorm(n =I,mu=rep(0,J),Sigma = diag(1,J))
                                
                                # 2. U, D, V obtained through performing SVD
                                svd1  =  svd(X)
                                
                                # 3.
                                W = svd1$v[,1:R]
                                
                                # 4.  Replace some elements of W by 0
                                if (propsparse !=0){
                                  W = apply(W, 2, function(x){
                                    t = quantile(abs(x),propsparse)
                                    x[abs(x) < t] = 0
                                    return(x)
                                  })
                                }
                                
                                # 5. Normalize each columns of W to a unit vector
                                W  <- W %*% diag(1/sqrt(diag(t(W) %*% W)))
                                
                                # 6. 
                                Z =  X%*%W # Component scores
                                #P = t(X)%*%ginv(t(Z)) 
                                Psvd = svd((t(X)%*%X)%*%W)
                                P = Psvd$u%*%t(Psvd$v)
                                # 7. 
                                Xtrue =  Z%*%t(P)
                                
                                
                                # Adding noice
                                SSqXtrue =  sum(Xtrue^2)                        # sum squares of the data set
                                EX = mvrnorm(I,mu= rep(0,J),Sigma = diag(1,J))  # EX = Error of X
                                SSqEX = sum(EX^2)                               # Sum squares fo the EX
                                fx = sqrt(SSqXtrue*(1-vafx)/(vafx * SSqEX))
                                Xnew = Xtrue + fx*EX                            # Data with noise
                                
                                # Saving data
                                out = list(X = Xnew, W =W ,P = P, Z = Z, k = R, 
                                           Propsparse = propsparse, nzeros =  n_zeros, teller = i)
                                save(out, file=paste0("DATA-R/Wsparse",i, ".RData"))
                              }
# stop Cluster
stopCluster(c1)
