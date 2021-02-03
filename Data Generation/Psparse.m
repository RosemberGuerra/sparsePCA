% P sparse generation using parfor %
% Rosember Guerra %
% 15-09-2019 %

clear
sc = parallel.pool.Constant(RandStream('Threefry'));

mkdir DATA-Matlab

% Define factors and levels
VAFx = [.80 .95 1];     % Variance accounted for X
p_sparse = [0 .5 .8];   % proportion of sparsity in the data
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
    
    X = randn(I,J);
    [T, s , P]=svds(X,R);               % Orthogonal T
    if propsparse ~= 0
        v = quantile(abs(P),propsparse);
        P(abs(P) < v) = 0;
    end
    P = P*s;
    Xtrue = T*P';
    W = pinv(Xtrue)*T;  % W estimation
    
    % Adding noice
    SSqXtrue = sum(sum(Xtrue.^2));      % sum squares of the data set
    EX = mvnrnd(zeros(J,1), eye(J),I);  % EX = Error of X
    SSqEX = sum(sum(EX.^2));            % Sum of squares of EX
    fx =  sqrt(SSqXtrue*(1-vafx)/(vafx*SSqEX)); % Rescale noise to desired level
    Xnew = Xtrue + fx*EX; %
    
    % Saving the data
    teller = i;
    parsave(sprintf('DATA-Matlab/Psparse%d.mat',i),Xnew,W,P,T,R,...
        propsparse,vafx,teller);
end
%% Auxiliar functions %%

function parsave(fname,Xnew,W,P,T,R,propsparse,vafx,teller)
save(fname, 'Xnew','W','P','T','R', 'propsparse','vafx','teller');
end