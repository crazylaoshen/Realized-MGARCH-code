function [c,ceq] = rbekk_constraint_targeting(theta) %#ok<*INUSL>
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
k =4;
numparams = k*k;
B = reshape(theta(1:numparams),k,k);
A = reshape(theta(numparams+(1:numparams)),k,k) ;
ceq = [];

m = kron(A,A) + kron(B,B);
c = abs(eig(m)) - .99998;
end
