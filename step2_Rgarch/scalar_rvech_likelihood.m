function [ll,lls,Ht]=scalar_rvech_likelihood(parameters,Y,X)
% Log likelihood for SCALAR_VT_VECH(P,O,Q) estimation
%
% USAGE:
%   [LL, LLS, HT] = scalar_vt_vech_likelihood(PARAMETERS,DATA,DATAASYM,P,O,Q,C,CASYM,KAPPA,BACKCAST,BACKCASTASYM,ISJOINT,USECOMPOSITE,ESTIMFLAG)
%
% INPUTS:
%   PARAMETERS   - A vector of vech GARCH process parameters: [alpha beta]' or
%                    [vech(C)' alpha beta]' if ISJOINT
%   DATA         - K by K by T matrix of covariance innovations
%   DATAASYM     - K by K by T matrix of asymmetric covariance innovations
%   P            - Positive, scalar integer representing the number of lags of the innovation process
%   O            - Non-negative scalar integer representing the number of lags of asymmetric process
%   Q            - Non-negative scalar integer representing the number of lags of conditional covariance
%   C            - The unconditional covariance of the data (cov(data)
%   CASYM        - The unconditional expectation of the asymmetric covariance
%   BACKCAST     - Back cast value for starting the recursion
%   BACKCASTASYM - Back cast value (asymetric terms) for starting the recursion
%   ISJOINT      - Logical indicating whether the likelihood is the joint across all parameters
%                    (intercept and dynamics) or is just for the dynamics
%   USECOMPOSITE - Indicates whether to use composite likelihood:
%                    0 - Use standard likelihood (not composite)
%                    1 - Use diagonal composite likelihood
%                    2 - Use full composite likelihood
%   ESTIMFLAG    - [OPTIONAL] Flag (0 or 1) to indicate if the function
%                    is being used in estimation.  If it is 1, then the parameters are transformed
%                    from unconstrained values to constrained by standard garch model constraints
%
% OUTPUTS:
%   LL           - Minus 1 times the log likelihood
%   LLS          - Time series of log likelihoods (Also multiplied by -1)
%   HT           - Time series of conditional covariances
%
% COMMENTS:
%   See also SCALAR_VT_VECH

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 4    Date: 10/28/2009

[k,~,T]=size(X);
%If nargin~=7, then we don't need to transform the parameters


beta=parameters(1);
alpha=parameters(2);
%Compute the constant

M = mean(X,3); % unconditional means of Realized kernel
% k = P\M;  % used in targeting parameter
C=cov(Y);% unconditional covariance matrix is just outproduct of returns
const=(C-alpha*M-beta*C);
%Initialize the covariance
Ht=zeros(k,k,T);

%Initialize the log likelihood
lls=zeros(T,1);

%Perform the recursion
 Ht(:,:,1)=const;
for t=2:T;
   
    Ht(:,:,t)=const+beta*Ht(:,:,t-1)+alpha*X(:,:,t-1);
    
end
% Replace these lines to make it work better with poorly conditioned covairance
% likelihoods(i)=likconst+(log(det(Ht(:,:,t)))+data(i,:)*Ht(:,:,t)^(-1)*data(i,:)');
% This is a trick to ensure better numerical stability
for t=1:T
    Q=sqrt(diag(Ht(:,:,t)));
    R=Ht(:,:,t)./(Q*Q');
    stdresid=Y(t,:)'*Y(t,:)./(Q*Q');
    lls(t)=0.5*(2*sum(log(Q))+log(det(R))+trace(R^(-1)*stdresid));
%     lls(t) = 0.5*log(det(Ht(:,:,t))) +0.5*sum(diag(Ht(:,:,t)\(Y(t,:)'*Y(t,:))));
% Wishart distribution
end
    
ll = sum(lls);