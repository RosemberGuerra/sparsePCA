PCAsimplimax <- function(X,k,ps){
  # X is the data matrix
  # ps is the proportion of sparsity of the loadings
  
  ## computing de loadings ##

  X = scale(X, center = TRUE, scale = FALSE) #

  Xpca =  svd(X, nv = k)
  Kloadings = Xpca$v%*%diag(Xpca$d[1:k])

  rloadings = as.matrix(Kloadings)
  Z = Xpca$u[,1:k]
  
  if (ps != 0){
    Nzeros =  ceiling(quantile(0:dim(X)[2],probs = ps))
    rotation = GPArotation::simplimax(Kloadings, k=Nzeros)
    rloadings = as.matrix(rotation$loadings)
    Z = Xpca$u[,1:k]%*%rotation$Th # rotation scores
    rloadings = apply(rloadings, 2, function(x){
      t = quantile(abs(x),probs = ps)
      x[abs(x) <= t] = 0
      return(x)
    })
    # rloadings =  apply(rloadings, 2, function(x){x/norm(x,"2")})
  }

  
  return(list(P = rloadings, Z = Z))
}