%% Application Gene expression data set for Autism %%
% Method: Gpower %
% Description: The index of sparseness (SI) is calculated and the 
% percentage of explained variance (PEV).
% Author: Rosember Guerra %
% Initial version: 03-05-2020 
% Edited version: 
%%
clear
close
clc

% load the data %
T = readtable('Data_Autism.csv', 'ReadRowNames',true);   
Autism= T{:,:};
X = Autism';  % The data is already processed, therefore, there is not need
% to scale it. 
clear Autism T
%% PCA %
[~,J] = size(X);    % Dimension of X.
[U,S,V] = svd(X);

% plot %
pev = diag(S).^2/sum(diag(S).^2);
PevK = cumsum(pev);

plot(PevK)

%% Plot 
c_gray = gray;
% colormap(gray)
% colormap(bone)
% colormap(colorcube)

c = [40*ones(1,7) 50*ones(1,6) 30*ones(1,14)]-10;
CG = c_gray(c,:);
scatter(U(:,1),U(:,2),[],c);
figure
scatter3(U(1:7,1),U(1:7,2),U(1:7,3),'*')%,CG,'filled')
hold on
scatter3(U(8:13,1),U(8:13,2),U(8:13,3),'+')%,CG,'filled')
scatter3(U(14:27,1),U(14:27,2),U(14:27,3),'o')%,CG,'filled')
xlabel('C. Score 1');
ylabel('C. Score 2');
zlabel('C. Score 3');
title('K=3, PCA')
%%
K = 3;          % Number of components.

upb = 0.76;     % This value was obtained manualy to get 
% 97.52% of proportion of sparsity with 2 components. 
dW = J*K;
N = 1000;
l = linspace(0,upb,N);
  
% 2 components %
Xpca = U(:,1:K)*S(1:K,1:K)*V(:,1:K)';
Rv0 = 1 - (norm(X-Xpca,'fro')/norm(X,'fro'))^2;

% initial values %
IS = zeros(N,1);
Rva = zeros(N,1);
ProS = zeros(N,1);
Upca = U(:,1:K);
%%
parfor i= 1:N
    %disp(i)
    gamma = l(i)*ones(1,K);
    W = GPower2(X,gamma,K,'l1',0,Upca);
    PS =  sum(W(:) == 0)/dW;
    
    T = X*W;
    P = X'*pinv(T');
    Xs = T*P';
    Rva(i) = 1 - (norm(X-Xs,'fro')/norm(X,'fro'))^2;
    ProS(i) = PS;
    IS(i) = Rv0* Rva(i)*PS;
    
end
%%
[~,Indmx] = max(IS);

hold off
plot(ProS,IS)
plot(ProS,Rva)
Results = table(ProS,IS,Rva,...
    'VariableNames',{'PS_Gpower' ,'IS_Gpower','PEV_Gpower'});
writetable(Results,'Result-Gpower_3.txt','Delimiter',',')
%%% Saving the vales %%%

gamma = l(Indmx)*ones(1,K);
W = GPower2(X,gamma,K,'l1',0,Upca);
PS =  sum(W(:) == 0)/dW;

T = X*W;
P = X'*pinv(T');
Xs = T*P';
Rvafinal = 1 - (norm(X-Xs,'fro')/norm(X,'fro'))^2;
save('AutismApplication_3.mat','W','T','U','S','X')

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

c = [40*ones(1,7) 50*ones(1,6) 30*ones(1,14)]-10;
CG = c_gray(c,:);
scatter(T(:,3),T(:,2),[],c);