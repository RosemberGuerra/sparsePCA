RelativeError = function(A,B){
  #  This function calculate the relative error between the columns of the
  #  matrix A and B.
  #  Rosember Guerra 
  #  12-09-2019 
  ##  Inputs
  #  A is the population matrix
  #  B is the matrix to be compared
  ##  Output
  #  er relative error
  #  p is the permutation that gives the max tc
  #  s is the sign of the 
 
  # Number of components
  k=dim(A)[2] 
  
  # Sign combination matrix
  Sign = c(1,-1)
  if (k==2){
    Sign_matrix = expand.grid(Sign,Sign)
  }else{
    Sign_matrix = expand.grid(Sign,Sign,Sign)
  }
  
  # Columns permutations
  Perms = combinat::permn(1:k)
  NPerms = length(Perms)
  Error_sign = matrix(rep(0,2*NPerms),ncol = 2)
  
  for (i in 1:NPerms) {
    PB = B[,Perms[[i]]]
    EP = rep(0,dim(Sign_matrix)[1])
    for (j in 1:dim(Sign_matrix)[1]) {
      PBS = PB%*%diag(Sign_matrix[j,])
      EP[j] = (norm((A-PBS),type = 'F')/norm(A,type = "F"))^2
    }
    Error_sign[i,1] = min(EP)
    Error_sign[i,2] = which(EP == Error_sign[i,1])[1]
  }
  R_error = min(Error_sign[,1])
  R_error_indx = which(Error_sign[,1] == R_error)[1]
  p = Perms[[R_error_indx]]
  s = Sign_matrix[Error_sign[R_error_indx,2],]
  return(list(Error = R_error, permutation = p, s= s))
}