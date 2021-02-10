%% Application Gene expression data set for Autism %%
% Method: sPCA-rSVD %
% Description: The index of sparseness (SI) is calculated and the 
% percentage of explained variance (PEV).
% Author: Rosember Guerra %
% created: 16-11-2020 
% Edited: 

% This application was implemented in a supercomputer with the following 
% characteristics:
% OS version: Windows 2012 R2
% CPU: 24 Core 2.60 GHz Intel Xeon(R)
% Memory: 523780 MB
%%
clear
close
clc

% load the data %
T = readtable('Data_Autism.csv', 'ReadRowNames',true);   
Autism= T{:,:};
X = Autism';  % The data is already processed, therefore, there is not need
[I,J]= size(X);
% to scale it. 
clear Autism T

K = 3;          % Number of components.
dP = J*K;
N = 1000;
l = linspace(0,0.98,N);
  
% 2 components %
%% PCA %
[~,J] = size(X);    % Dimension of X.
[U,S,V] = svd(X);
% [Ui,Si,Vi] = svd(X);
Uini = U(:,1:K);
Vini = V(:,1:K)*S(1:K,1:K);
Sini= S(1:K,1:K);

Xpca = U(:,1:K)*S(1:K,1:K)*V(:,1:K)';
Rv0 = 1 - (norm(X-Xpca,'fro')/norm(X,'fro'))^2;

%% sPCArSVD 

[P,T] = sPCArSVD2(X,K,.943,Uini,Vini,Sini);

PS =  sum(P(:) == 0)/dP;

Xs = T*P';
Rvafinal = 1 - (norm(X-Xs,'fro')/norm(X,'fro'))^2;
% save('AutismApplication_spcarsvd.mat','P','T','U','S','X')

% plot %

scatter3(T(1:7,1),T(1:7,2),T(1:7,3),'k*')%,CG,'filled')
hold on
scatter3(T(8:13,1),T(8:13,2),T(8:13,3),'kd','filled')%,CG,'filled')
scatter3(T(14:27,1),T(14:27,2),T(14:27,3),'ko')%,CG,'filled'
legend({'dup15','FMR1','Control'},'Location','northwest')
xlabel('C. Score 1');
ylabel('C. Score 2');
zlabel('C. Score 3');

hold off
