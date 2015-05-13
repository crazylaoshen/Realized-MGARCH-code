function [c,ceq] = constraint_diagonal(theta) %#ok<*INUSL>
% Non-linear constraint for estimation of BEKK(p,o,q) multivariate volatility models
%
% USAGE:
%  [C,CEQ] = bekk_constraint(PARAMETERS,DATA,DATAASYM,P,O,Q,BACKCAST,BACKCASTASYM,TYPE) 
%
% INPUTS:
%   See bekk_likelihood
%
% OUTPUTS:
%   C   - Vector of inequality constraints
%   CEQ - Empty
%
% COMMENTS:
%
%  EXAMPLES:
%
% See also BEKK, BEKK_LIKELIHOOD

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 1    Date: 3/27/2012
% k = 2;
k =10;
numParams = k;
k2 = k*(k+1)/2;
tempP = theta(k2+(1:numParams));
B = diag(tempP);
% k2 = k2+numParams;
% tempP = theta(k2+(1:numParams));
% A = reshape(tempP,k,k);
ceq = [];
m = kron(B,B);
% m = kron(A,A) + kron(B,B);
c = abs(eig(m)) - .99998;
end
