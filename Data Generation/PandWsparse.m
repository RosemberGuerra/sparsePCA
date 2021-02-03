% P and W sparse data generation using parfor
% Rosember Guerra
% 15-09-2019

clear
sc = parallel.pool.Constant(RandStream('Threefry'));

mkdir DATA-Matlab

% Define factors and levels
VAFx = [.80 .95 1];     % Variance accounted for X
p_sparse = [0.7 .8 .9];   % proportion of sparsity in the data
s_size = [100 500];     % Sample size
n_variables = [10 100 1000];% Number of variables
n_components = [2 3];   % Number of components
n_replications = 100;   % Number of replications

% design matrix

% desig_matrix = setting_matrix(n_variables,s_size,p_sparse,n_components,VAFx);
[grid_Nvariables,grid_Ssize, grid_psparse,   grid_Ncomponents, grid_VAFx]=...
    ndgrid(n_variables,s_size,p_sparse,n_components,VAFx);

desig_matrix = [grid_Nvariables(:) grid_Ssize(:) grid_psparse(:)...
   grid_Ncomponents(:)  grid_VAFx(:)];

desig_matrix_repetitions = repmat(desig_matrix,n_replications,1);

Name = {'n_variables','s_size','p_sparse','n_components','VAFx' };

Table = table(desig_matrix_repetitions(:,1),desig_matrix_repetitions(:,2),...
    desig_matrix_repetitions(:,3),desig_matrix_repetitions(:,4),...
    desig_matrix_repetitions(:,5),'VariableNames',Name);

Ndatasets = size(desig_matrix_repetitions,1);
save('DATA-Matlab/Info_simulation.mat','Ndatasets','desig_matrix_repetitions',...
    'desig_matrix','n_replications','Name','Table')

parfor i=1:size(desig_matrix_repetitions,1)
    % Random generation %
    stream = sc.Value;        % Extract the stream from the Constant
    stream.Substream = i;
    % set the specific values of the parameters
    variables = num2cell(desig_matrix_repetitions(i,:));
    [J, I, propsparse, R,vafx] = deal(variables{:});
    
    %
    nonzeros = floor(quantile(1:J,1-propsparse));
    X = randn(I,J);

    [U, s , V]=svds(X,R);  
    
    % 7.
    indxnonzero =  reshape(randsample(J,nonzeros*R),[nonzeros R]);
    zeroindx = ones(J,R);
    for z=1:R
        zeroindx(indxnonzero(:,z),z) = 0;
        V(zeroindx(:,z)==1,z) = 0;
    end
    % P and W are equal under this specification %
    V = V*diag(1./sqrt(diag(V'*V)));
    W2 = V; 
    %% 
    P = V*s; % This P2 is used for sPCA-rSVD method, for fair comparison.
    T = U;
    T2 = X*W2;
%     P2 = X'*pinv(T2');
    %%
    Xtrue = T*P';
    
    Xtrue2 = T2*W2';
    % Adding noice 1 %
    
    SSqXtrue = sum(sum(Xtrue.^2));      % sum squares of the data set
    EX = mvnrnd(zeros(J,1), eye(J),I);  % EX = Error of X
    SSqEX = sum(sum(EX.^2));            % Sum of squares of EX
    fx =  sqrt(SSqXtrue*(1-vafx)/(vafx*SSqEX)); % Rescale noise to desired level
    
    Xnew = Xtrue + fx*EX; % adding additive noise!
    
    % Adding noice 2 %
    
    SSqXtrue = sum(sum(Xtrue2.^2));      % sum squares of the data set
    EX = mvnrnd(zeros(J,1), eye(J),I);  % EX = Error of X
    SSqEX = sum(sum(EX.^2));            % Sum of squares of EX
    fx =  sqrt(SSqXtrue*(1-vafx)/(vafx*SSqEX)); % Rescale noise to desired level
    
    Xnew2 = Xtrue2 + fx*EX; % adding additive noise!
    
    % Saving the data 
    teller = i;
    parsave(sprintf('DATA-Matlab/PWsparse%d.mat',i),Xnew,Xnew2,P,W2,T,T2,R,...
        propsparse,vafx,teller);
    
    
end
%% Auxiliar functions %%

function parsave(fname,Xnew,Xnew2,P,W2,T,T2,R,propsparse,vafx,teller)
save(fname, 'Xnew','Xnew2','P','W2','T','T2','R', 'propsparse','vafx','teller');
end