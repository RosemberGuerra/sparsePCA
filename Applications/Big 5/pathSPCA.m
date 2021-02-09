function W = pathSPCA(X,k,ps)
% S should be the covariance matrix         %
% k is the number components                %
% ps proportion of sparsity                 %
% Ik is the index with no zero variables    %
% W weighst matrix                          %

[~,p]=size(X); 
% S = X'*X/n;             %Covariance matrix
S = cov(X);


if ps == 0 % 0% sparsity
    [V,D] = eig(S);
    [~,ind] = sort(diag(D),'descend');
    
    W = V(:,ind(1:k));
    
else 
    NonZeros = floor(quantile(1:p,1-ps)); % cardinality
    S = PD(S);
    [~,ix] = sort(diag(S),'descend');
    S = S(ix,ix);                   % Sort covariance w.r.t. diag
    [~,~,res1]=FullPathCov(S);
    ik = res1(1:NonZeros,NonZeros); % Nonzero index
    S = S(ik,ik);
    [V,D] = eig(S);
    [~,ind] = sort(diag(D),'descend');
    V = V(:,ind(1));
    W = zeros(p,1);
    W(ik) = V;
    W = W(ix);  % Back to the original order !!!

end

if k > 1  % More than one component
    
    Xres = X- X*(W*W')/norm(W)^2;
    kres = k-1;
    Wres = pathSPCA(Xres,kres,ps); % Recursive function
    W = [W Wres];
end

end