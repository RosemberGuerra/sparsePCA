function rv = correlation(A,B)
% This function calculate the Tucker Congruence
% Rosember Guerra %
% 13-09-2019 %
%% Inputs
% A is the population matrix
% B is the matrix to be compared
% A and B should have the same dimension
%% Output
% rv is the Tucker congruence
%%
k = size(A,2);
Numerator = abs(sum(A.*B));
Denominator = sqrt(sum(A.*A).*sum(B.*B));
rv = sum(Numerator./Denominator)/k;
end