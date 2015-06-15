function [c,ceq] = constraint_rotate(theta) %#ok<*INUSL>
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

tempP = theta(1:k);
B = diag(tempP);
tempP = theta(k+(1:k));
A = diag(tempP);
ceq = [];
% c = diag((A.^2) + (B.^2))-1;
m = kron(A,A) + kron(B,B);
c = eig(m)-0.99998;

end
