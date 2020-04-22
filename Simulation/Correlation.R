correlation = function(A,B){
  
  k = dim(A)[2]
  Numerator = abs(colSums(A*B));
  Denominator = sqrt(colSums(A*A)*colSums(B*B))
  rv = sum(Numerator/Denominator)/k;
  
}