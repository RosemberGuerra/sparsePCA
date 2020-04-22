function [er,p,s] = ErrorRelative(A,B)
% This function calculate the Tucker Congruence
% Rosember Guerra %
% 12-09-2019 %
%% Inputs
% A is the population matrix
% B is the matrix to be compared
%% Output
% er relative error
% p is the permutation that gives the max tc
% s is the sign of the combiantion
%%
k = size(A,2); % Number of components
% k has only to dimension 2 and 3.
% Possible sign
if k==2
    [S1,S2] = ndgrid([1 -1],[1 -1]);
    Sign_matrix = [S2(:) S1(:)];
else 
    [S1,S2,S3] = ndgrid([1 -1],[1 -1],[1 -1]);
    Sign_matrix = [S3(:) S2(:) S1(:)]; 
end

Perms = perms(1:k); % all possible combination
NPerms = size(Perms,1);
Error_sign = zeros(NPerms,2);

for i=1:NPerms
    PB = B(:,Perms(i,:));
    EP = zeros(size(Sign_matrix,1),1);
    for j=1:size(Sign_matrix,1)
        PBS = PB*diag(Sign_matrix(j,:));
        EP(j) = (norm(A-PBS,'fro')/norm(A,'fro'))^2;
    end
    [Error_sign(i,1), Error_sign(i,2)] = min(EP);
    
end

[er, MinIdx] = min(Error_sign(:,1));
p = Perms(MinIdx,:);
s = Sign_matrix(Error_sign(MinIdx,2),:);


end