%% Simulation pathSPCA parfor%%
% Rosember Guerra %
% Initial version: 09-10-2019  
% Edited version: 


clear
close
clc

% simulation %
% 1 means the model that match the data set     %
% 2 means P and W sparse and equals             %
% 3 means the model and the data don't match    %

% setting directory-W sparse
DataDir1 = '../../DataGeneration/Wsparse/DATA-Matlab/';

% Setting directories- P and W sparse %
DataDir2 = '../../DataGeneration/PandWsparse/DATA-Matlab/';

% Setting directories- Cross data %
DataDir3 = '../../DataGeneration/Psparse/DATA-Matlab/';

load([DataDir1,'Info_simulation.mat'])
%% Performance measure variables %%
% 1.
MR01 = zeros(Ndatasets,1);           
ErrorLW1 = zeros(Ndatasets,1);  
ErrorScores1 = zeros(Ndatasets,1);   

Svar1  =  zeros(Ndatasets,1);           % Sparse variance
% 2.
MR02 = zeros(Ndatasets,1);           
ErrorLW2 = zeros(Ndatasets,1);  
ErrorScores2 = zeros(Ndatasets,1);   

Svar2  =  zeros(Ndatasets,1);           % Sparse variance

% 3.
MR03 = zeros(Ndatasets,1);           
Correlation3 = zeros(Ndatasets,1);
CorrelationZ3 = zeros(Ndatasets,1);

Svar3  =  zeros(Ndatasets,1);           % Sparse variance

%% -------------------------------- %%
parfor i=1:Ndatasets
    %% Estimation %%
    % Loading data
    Out1 = load([DataDir1,sprintf('Wsparse%d.mat',i)]);     % W sparse
    Out2 = load([DataDir2,sprintf('PWsparse%d.mat',i)]);    % P and W sparse
    Out3 = load([DataDir3,sprintf('Psparse%d.mat',i)]);     % P sparse
    
    
    W1 = pathSPCA(Out1.Xnew,Out1.R,Out1.propsparse);
    W2 = pathSPCA(Out2.Xnew2,Out2.R,Out2.propsparse);
    W3 = pathSPCA(Out3.Xnew,Out3.R,Out3.propsparse);
    % Component scores %
    T1 = Out1.Xnew*W1; 
    T2 = Out2.Xnew2*W2;
    T3 = Out3.Xnew*W3;

    %% Performance Measurements %%
    % Relative error %
    [ErrorLW1(i),perm1,s1] = ErrorRelative(Out1.W ,W1 );
    [ErrorLW2(i),perm2,s2] = ErrorRelative(Out2.W2 ,W2 );
    [~,perm3,s3] = ErrorRelative(Out3.W ,W3 );
    %
    W1 = W1(:,perm1)*diag(s1);
    W2 = W2(:,perm2)*diag(s2);
    W3 = W3(:,perm3)*diag(s3);
    %
    T1 = T1(:,perm1)*diag(s1);    
    T2 = T2(:,perm2)*diag(s2);    
    T3 = T3(:,perm3)*diag(s3);
    %
    ErrorScores1(i) = (norm(Out1.T-T1,'fro')/norm(Out1.T,'fro'))^2;
    ErrorScores2(i) = (norm(Out2.T2-T2,'fro')/norm(Out2.T2,'fro'))^2;
    Correlation3(i) = correlation(Out3.P, W3);
    CorrelationZ3(i) = correlation(Out3.T, T3);
    
    % Misidentificade rate %
    ZeroInd1 = find(Out1.W(:)== 0); % W sparse
    ZeroInd2 = find(Out2.W2(:)== 0); % P and W sparse
    ZeroInd3 = find(Out3.P(:)== 0); % P sparse
    if Out1.propsparse == 0
        MR01(i) = 0;
        MR03(i) = 0;
    else
        wz1 = W1(ZeroInd1);
        MR01(i) =1- sum(wz1==0)/length(ZeroInd1);
        wz3 = W3(ZeroInd3);
        MR03(i) =1- sum(wz3==0)/length(ZeroInd3);
    end
    if Out2.propsparse == 0
        MR02(i) = 0;        
    else
        wz2 = W2(ZeroInd2);
        MR02(i) =1- sum(wz2==0)/length(ZeroInd2);
    end
        
    
    
    % Variance explained-2 %
    P1 = Out1.Xnew'*pinv(T1');
    Xs1 = T1*P1';
    Svar1(i) = 1 - (norm(Out1.Xnew-Xs1,'fro')/norm(Out1.Xnew,'fro'))^2;
    P2 = Out2.Xnew2'*pinv(T2');
    Xs2 = T2*P2';
    Svar2(i) = 1 - (norm(Out2.Xnew2-Xs2,'fro')/norm(Out2.Xnew2,'fro'))^2;
    P3 = Out3.Xnew'*pinv(T3');
    Xs3 = T3*P3';
    Svar3(i) = 1 - (norm(Out3.Xnew-Xs3,'fro')/norm(Out3.Xnew,'fro'))^2;
    
    
end
%% Joining the data set %%
Methodology = repmat('pathSPCA',Ndatasets, 1);
PW_info = load([DataDir2,'Info_simulation.mat']);
p_sparse2 = PW_info.desig_matrix_repetitions(:,3);

Table2 = table(ErrorLW1,ErrorScores1,MR01,Svar1,...
    ErrorLW2,ErrorScores2,MR02,Svar2,...
    Correlation3,CorrelationZ3,MR03,Svar3);

DataPathSPCA = [table(Methodology,p_sparse2) Table Table2];
writetable(DataPathSPCA,'JointDataPathSPCA.txt','Delimiter',',') 
save('ResultsPathSPCA.mat')
%% Auxiliar functions %%
% function parsave(fname,W,T)
%     save(fname,'W','T');
% end