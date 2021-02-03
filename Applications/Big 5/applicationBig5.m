%% Application Big5 data set %%
% Method: sPCArSVD %
% Description: The index of sparseness (SI) is calculated and the 
% percentage of explained variance (PEV).
% Author: Rosember Guerra %
% Initial version: 21-01-2020 
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

plot(PevK)

% 
nz = length(NZeros);
IS = zeros(nz,1);
Rva = zeros(nz,1);

% 
for i= 1:nz
% for i= 60:nz
    disp(i)
    
    %
    [P,T] = sPCArSVD(X,K,PS(i));
    Xs = T*P';
    Rva(i) = 1 - (norm(X-Xs,'fro')/norm(X,'fro'))^2;
    IS(i) = Rv0* Rva(i)*PS(i);
    
end
plot(PS,IS)



%%% Saving the vales %%%
Results = table(PS',IS,Rva,...
    'VariableNames',{'PS' ,'IS_spcarsvd','PEV_spcarsvd'});
writetable(Results,'Result-sPCArSVD.txt','Delimiter',',') 

% Taking the max IS %
[~,Indx] = max(IS);

[P,T] = sPCArSVD(X,K,PS(Indx));
Xs = T*P';


save('sPCArSVDApplication.mat','X','Xs','T','P')

