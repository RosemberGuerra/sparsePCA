%% Application Big5 data set for resubmission %%
% Method: GPower and pathspca %
% Description: The index of sparseness (SI) is calculated and the 
%   percentage of explained variance (PEV).
% Author: Rosember Guerra-Urzola%
% Initial version: 16-11-2020 
% Edited version: 
%%
clear
close
clc

load('big5.txt') % Loading the data set.


X = normalize(big5); %BigFive;

[~,J] = size(X);    % Number of variables

K = 5;          % Number of components.
% ncomp = length(K);

NZeros =  0:(J-1); 
PS = NZeros/J;

% PCA %
[U,S,V] = svd(X);

% components %
Xpca = U(:,1:K)*S(1:K,1:K)*V(:,1:K)';
Rv0 = 1 - (norm(X-Xpca,'fro')/norm(X,'fro'))^2;

% diagS = diag(S);
% sum(diagS(1:K).^2)/sum(diagS.^2)

% plot %
pev = diag(S).^2/sum(diag(S).^2);
PevK = cumsum(pev);

Wgp=GP(X,K,.734,1);

Tgp=X*Wgp;
Pgp = X'*pinv(Tgp');
Xs_gp = Tgp*Pgp';
Rva_gp = 1 - (norm(X-Xs_gp,'fro')/norm(X,'fro'))^2;

Wpth=pathSPCA(X,K,.734);
Tpth=X*Wpth;
Ppth = X'*pinv(Tpth');
Xs_pth = Tpth*Ppth';
Rva_pth = 1 - (norm(X-Xs_pth,'fro')/norm(X,'fro'))^2;


ocean_gp= zeros(5,5);
ocean_pth= zeros(5,5);
N=1:5:240;


% NEOAC order % 35241

for i=1:5
    for j=1:5
        ind= N+j-1;
        ocean_gp(i,j)=sum(Wgp(ind,i)~=0);
        ocean_pth(i,j)=sum(Wpth(ind,i)~=0);
    end
end

ocean_gp =ocean_gp(:,[3 5 2 4 1]);
ocean_pth =ocean_pth(:,[3 5 2 4 1]);