function [V,U]= sPCArSVD(X,k,ps)
%% sPCA-rSVD %%
% Rosember Guerra 
% 20-05-2019 
% Description: This function calculate the component scores and component 
% loadings using Shen and Yuan (2008) methodology. For the tuning-parameter 
% is adjusted for getting the proportion of sparsity as described pag. 1020 on the paper. 

%% Inputs:
% X: data matrix, it must be centered matrix.
% k: Number of components.
% ps: Proportion of sparsity
% V0 initial value for V
% U0 initial value for U
%% Output:
% V: loadings vectors
% u: Component vectors
%%
rng('default'); % For random number generation
[~, p] = size(X); % Centering of the data
eps = 1e-7; % 
MaxIter = 1000;


dif = 1;
iter = 0;
% 1. SVD
[Ui,Si,Vi] = svd(X);
uold = Ui(:,1);
vold = Si(1,1)*Vi(:,1);
Nzeros = ceil(p*ps); % number of zeros
if ps == 0
    l = 0;
else
    AXu = sort(abs(X'*uold));
    l = AXu(Nzeros);
end
fold = norm(X-uold*vold')+2*l*sum(abs(vold));

% 2.  Update
while dif > eps && iter < MaxIter 
    iter = iter +1;
    % a.
    Xu = X'*uold;
    vnew = sign(Xu).*max(0,abs(Xu)-l);
    if sum(vnew == 0)~= Nzeros && ps ~= 0
        AXu = sort(abs(Xu));
        l = AXu(Nzeros);
        vnew = sign(Xu).*max(0,abs(Xu)-l);
    end
    
    % b.
    Xv = X*vnew;
    unew = Xv/norm(Xv);

    fnew = norm(X-unew*vnew')+2*l*sum(abs(vnew));
    dif = abs(fold-fnew)/(1+fold);
    uold = unew;
    fold = fnew; 
end

vnew = Si(1,1)*vnew/norm(vnew);

V1 =  vnew;
U1 =  unew;
if k > 1 % For more than 1 component
    Xres = X-unew*vnew'; % Residual matrix
    kres = k-1;
    [Vres,Ures]= sPCArSVD(Xres,kres,ps); % Recursive function
    V = [V1 Vres];
    U = [U1 Ures];
else
    V = V1;
    U = U1;
end

end