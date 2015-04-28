function [L, H] = rbekk_likelihood_scalar(theta, Y, X)
% code revised based on Dr. Denis Pelletier's course in asset pricing @ NCSU

% Minus the log-likelihood of a DCC model

% input: theta --- parameter composed by one intercept and two coefficient matrix
% each one has k(k+1)/2 theta, so theta has k(k+1) theta
 %       V  ---   matrix of returns
 %       X -----  matrix of realized mesure ,as RK

% k = size(Y,2);
% k2 = k*(k+1)/2;
% vechC = theta(1:k2);
% C=zeros(k);
% Use a logical trick to inverse the vech
% ind=tril(true(k));
% C(ind)=vechC;
% C = C*C';
B = theta(1);
A = theta(2);

 
P = cov(Y);% Y is lret, P is COVARIANCE MATRIX FOR STANDARDIZED RESIDUALS
P = logm(P);
M = mean(X,3); % unconditional means of Realized kernel
% k = P\M;  % used in targeting parameter

[K,~,T] = size(X);
logH = zeros(K,K,T);
H = zeros(K,K,T);
logH(:,:,1)=P; % we set the unconditional covariance matrix of daily(close to close) return as starting point of revolved conditional covariance
log_likelihood = zeros(T,1);
intercept = (P- B*P -A*M);

for t=2:T
  logH(:,:,t) = intercept + B*logH(:,:,t-1)+ A*X(:,:,t-1) ; %conditional covariance MATRIX, insure positive semidefinite
end

% Compute minus log-likelihood
for t=1:T
        H(:,:,t) = expm(logH(:,:,t));
    Q=sqrt(diag(H(:,:,t)));
    R=H(:,:,t)./(Q*Q');
    %stdresid=X(:,:,t)./(Q*Q');
    stdresid=Y(t,:)'*Y(t,:)./(Q*Q');
    log_likelihood(t)=0.5*(2*sum(log(Q))+log(det(R))+trace(R^(-1)*stdresid));
% lls(t) = 0.5*log(det(H(:,:,t))) +0.5*sum(diag(H(:,:,t)\(Y(t,:)'*Y(t,:))));% Wishart distribution
%   log_likelihood(t) = 0.5*log(det(H(:,:,t))) +0.5*trace((H(:,:,t))^(-1)*(Y(t,:)'*Y(t,:)));% Wishart distribution
end

L = sum(log_likelihood);