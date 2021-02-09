function Z = GP(X,k,ps,L)
%%----inputs----%%
% X: data matrix
% k: Number of components
% ps: proportion of sparsity
% L: pelanty value, should be a either 0 or 1
%----Outputs-----%%
% Z: Matrix of component weights
%%
[~,p] = size(X);
% X = X - repmat((mean(X,1)),n,1); % Centering data matrix
Nzeros = floor(p*k*ps); % number of zeros
nzeros = 0; % initial value for the number of zeros
MaxIter = 1000;

if L == 0 % Type of penalty
    l = 'l0';
    up = max(vecnorm(X));
else
    l = 'l1';
    up = max(vecnorm(X))^2;
end
if Nzeros == 0 % No sparsity 
    Z=GPower(X,zeros(1,k),k,l,0);
end
low = 0; % Initial values for binary search
iter = 0;

 while Nzeros ~= nzeros && MaxIter > iter
     iter = iter + 1;
     g = (low+up)/2;
     if L == 0
         gamma=(g^2)*ones(1,k);
     else
         gamma=g*ones(1,k);
     end
     
     Z = GPower(X,gamma,k,l,0);
     nzeros = sum(Z(:) == 0);
     if Nzeros > nzeros
        low = g;
     elseif Nzeros < nzeros
        up = g;
     end
 end
 
end