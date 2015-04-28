function [L, H] = rbekk_likelihood_nontargeting(theta, Y, X)
% code revised based on Dr. Denis Pelletier's course in asset pricing @ NCSU

% Minus the log-likelihood of a DCC model

% input: theta --- parameter composed by one intercept and two coefficient matrix
% each one has k(k+1)/2 theta, so theta has k(k+1) theta
 %       V  ---   matrix of returns
 %       X -----  matrix of realized mesure ,as RK

k = size(Y,2);
numParams = k*k;
k2 = k*(k+1)/2;
vechC = theta(1:k2);
C=zeros(k);
% Use a logical trick to inverse the vech
ind=tril(true(k));
C(ind)=vechC;
C = C*C';
tempP = theta(k2+(1:numParams));
B = reshape(tempP,k,k);
k2 = k2+numParams;
tempP = theta(k2+(1:numParams));
A = reshape(tempP,k,k);
 
P = cov(Y);% Y is lret, P is COVARIANCE MATRIX FOR STANDARDIZED RESIDUALS
% M = mean(X,3); % unconditional means of Realized kernel
% k = P\M;  % used in targeting parameter

[N,~,T] = size(X);
H = zeros(N,N,T);
H(:,:,1)=P; % we set the unconditional covariance matrix of daily(close to close) return as starting point of revolved conditional covariance
log_likelihood = zeros(T,1);

for t=2:T
  H(:,:,t) = C + B*H(:,:,t-1)*B'+ A*X(:,:,t-1)*A' ; %conditional covariance MATRIX, insure positive semidefinite
end

% Compute minus log-likelihood
for t=1:T
  log_likelihood(t) = 0.5*log(det(H(:,:,t))) +0.5*trace((H(:,:,t))^(-1)*(Y(t,:)'*Y(t,:)));% Wishart distribution
end

L = sum(log_likelihood);