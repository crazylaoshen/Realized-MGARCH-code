function [L, H] = rbekk_likelihood_diagonal(theta, Y, X)
% code revised based on Dr. Denis Pelletier's course in asset pricing @ NCSU

% Minus the log-likelihood of a DCC model

% input: theta --- parameter composed by one intercept and two coefficient matrix
% each one has k(k+1)/2 theta, so theta has k(k+1) theta
 %       V  ---   matrix of returns
 %       X -----  matrix of realized mesure ,as RK

[K,~,T] = size(X);
numParams = K;
k2 = K*(K+1)/2;
vechC = theta(1:k2);
C=zeros(K);
% Use a logical trick to inverse the vech
ind=tril(true(K));
C(ind)=vechC;
C = C*C';
tempP = theta(k2+(1:numParams));
B = diag(tempP);
k2 = k2+numParams;
tempP = theta(k2+(1:numParams));
A = diag(tempP);
 
P = cov(Y);% Y is lret, P is COVARIANCE MATRIX FOR STANDARDIZED RESIDUALS
P = logm(P);
% M = mean(X,3); % unconditional means of Realized kernel
% k = P\M;  % used in targeting parameter


H = zeros(K,K,T);
logH=zeros(K,K,T);
logH(:,:,1)=P; % we set the unconditional covariance matrix of daily(close to close) return as starting point of revolved conditional covariance
log_likelihood = zeros(T,1);

for t=2:T
  logH(:,:,t) = C + B*logH(:,:,t-1)*B'+ A*X(:,:,t-1)*A' ; %conditional covariance MATRIX, insure positive semidefinite
end

% Compute minus log-likelihood
parfor t=1:T
    H(:,:,t) = expm(logH(:,:,t));
    Q=sqrt(diag(H(:,:,t)));
    R=H(:,:,t)./(Q*Q');
    %stdresid=X(:,:,t)./(Q*Q');
    stdresid=Y(t,:)'*Y(t,:)./(Q*Q');
    log_likelihood(t)=0.5*(2*sum(log(Q))+log(det(R))+trace(R^(-1)*stdresid));
%   log_likelihood(t) = 0.5*log(det(H(:,:,t))) +0.5*trace((H(:,:,t))^(-1)*(Y(t,:)'*Y(t,:)));% Wishart distribution
end

L = sum(log_likelihood);