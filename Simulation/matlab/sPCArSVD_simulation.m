%% Simulation sPCArSVD parfor%%
% Rosember Guerra %
% Initial version: 18-09-2019 
% Edited version: 08-10-2019

clear
close
clc

% simulation %
% 1 means the model match the data set     	%
% 2 means P and W sparse and equals             %
% 3 means the model and the data don't match    %

% setting directory-P sparse
DataDir1 = '../../DataGeneration/Psparse/DATA-Matlab/';
% ResultsDir1 = 'Results/sPCArSVD1/';

% Setting directories- P and W sparse %
DataDir2 = '../../DataGeneration/PandWsparse/DATA-Matlab/';
% ResultsDir2 = 'Results/sPCArSVD2/';

% Setting directories- Cross data %
DataDir3 = '../../DataGeneration/Wsparse/DATA-Matlab/';
% ResultsDir3 = 'Results/sPCArSVD3/';

load([DataDir1,'Info_simulation.mat'])
%% Performance measure variables %%
% 1.
MR01 = zeros(Ndatasets,1);           
ErrorLW1 = zeros(Ndatasets,1);  
ErrorScores1 = zeros(Ndatasets,1);   

Svar1 = zeros(Ndatasets,1);
% 2.
MR02 = zeros(Ndatasets,1);           
ErrorLW2 = zeros(Ndatasets,1);  
ErrorScores2 = zeros(Ndatasets,1);   

Svar2 = zeros(Ndatasets,1);
% 3.
MR03 = zeros(Ndatasets,1);           
Correlation3 = zeros(Ndatasets,1);
CorrelationZ3 = zeros(Ndatasets,1);  

Svar3 = zeros(Ndatasets,1);
%%
parfor i=1:Ndatasets
    % Estimation %%
    % Loading data
    Out1 = load([DataDir1,sprintf('Psparse%d.mat',i)]);     % P sparse
    Out2 = load([DataDir2,sprintf('PWsparse%d.mat',i)]);    % P and W sparse
    Out3 = load([DataDir3,sprintf('Wsparse%d.mat',i)]);     % W sparse
    
    % 1.
%     DATA1 = Out1.Xnew;
%     [I1, ~] = size(DATA1); 
%     DATA1 = DATA1-repmat((mean(DATA1,1)),I1,1); 
    [P1,T1] = sPCArSVD(Out1.Xnew,Out1.R,Out1.propsparse);
    % 2.
%     DATA2 = Out2.Xnew;
%     [I2, ~] = size(DATA2); 
%     DATA2 = DATA2-repmat((mean(DATA2,1)),I2,1); 
    [P2,T2] = sPCArSVD(Out2.Xnew,Out2.R,Out2.propsparse);
    % 3.
%     DATA3 = Out3.Xnew;
%     [I3, ~] = size(DATA3); 
%     DATA3 = DATA3-repmat((mean(DATA3,1)),I3,1); 
    [P3,T3] = sPCArSVD(Out3.Xnew,Out3.R,Out3.propsparse);
    
    % Saving the results %
%     parsave(sprintf([ResultsDir1,'results%d.mat'],i),P1,T1);
%     parsave(sprintf([ResultsDir2,'results%d.mat'],i),P2,T2);
%     parsave(sprintf([ResultsDir3,'results%d.mat'],i),P3,T3);
%     
    % Performance Measurements %%
    % Relative error %
    [ErrorLW1(i),perm1,s1] = ErrorRelative(Out1.P ,P1 );
    [ErrorLW2(i),perm2,s2] = ErrorRelative(Out2.P ,P2 );
    [~,perm3,s3] = ErrorRelative(Out3.P ,P3 );
    
    P1 = P1(:,perm1)*diag(s1);
    P2 = P2(:,perm2)*diag(s2);
    P3 = P3(:,perm3)*diag(s3);
    
    T1 = T1(:,perm1)*diag(s1);
    T2 = T2(:,perm2)*diag(s2);
    T3 = T3(:,perm3)*diag(s3);
    
    ErrorScores1(i) = (norm(Out1.T-T1,'fro')/norm(Out1.T,'fro'))^2;
    ErrorScores2(i) = (norm(Out2.T-T2,'fro')/norm(Out2.T,'fro'))^2;
    %ErrorScores3(i) = (norm(Out3.T-T3,'fro')/norm(Out3.T,'fro'))^2;
    Correlation3(i) = correlation(Out3.W, P3);
    CorrelationZ3(i) = correlation(Out3.T, T3);
    
    % Misidentificade rate %
    ZeroInd1 = find(Out1.P(:)== 0); % P sparse
    ZeroInd2 = find(Out2.P(:)== 0); % P and W sparse
    ZeroInd3 = find(Out3.W(:)== 0); % W sparse
    if Out1.propsparse == 0
        MR01(i) = 0;
        MR03(i) = 0;
    else
        wz1 = P1(ZeroInd1);
        MR01(i) =1- sum(wz1==0)/length(ZeroInd1);
        wz3 = P3(ZeroInd3);
        MR03(i) =1- sum(wz3==0)/length(ZeroInd3);
    end
    if Out2.propsparse == 0
        MR02(i) = 0;        
    else
        wz2 = P2(ZeroInd2);
        MR02(i) =1- sum(wz2==0)/length(ZeroInd2);
    end
  
    
    
    
    % Sparse Explained variance 2 %
    Xs1 = T1*P1';
    Svar1(i) = 1 - (norm(Out1.Xnew-Xs1,'fro')/norm(Out1.Xnew,'fro'))^2;
    
    Xs2 = T2*P2';
    Svar2(i) = 1 - (norm(Out2.Xnew-Xs2,'fro')/norm(Out2.Xnew,'fro'))^2;
    
    Xs3 = T3*P3';
    Svar3(i) = 1 - (norm(Out3.Xnew-Xs3,'fro')/norm(Out3.Xnew,'fro'))^2;
end
%% Joining the data set %%
Methodology = repmat('sPCArSVD',Ndatasets, 1);
PW_info = load([DataDir2,'Info_simulation.mat']);
p_sparse2 = PW_info.desig_matrix_repetitions(:,3);


Table2 = table(ErrorLW1,ErrorScores1,MR01,Svar1,...
    ErrorLW2,ErrorScores2,MR02,Svar2,...
    Correlation3,CorrelationZ3,MR03,Svar3);


DatasPCArSVD = [table(Methodology,p_sparse2) Table Table2];
writetable(DatasPCArSVD,'JointDatasPCArSVD.txt','Delimiter',',') 
save('ResultsSPCArSVD.mat')

%% Auxiliar functions %%
% function parsave(fname,P,T)
% save(fname,'P','T');
% end