function [L, H] = rbekk_likelihood_targeting(theta, Y, X)
% code revised based on Dr. Denis Pelletier's course in asset pricing @ NCSU

% Minus the log-likelihood of a DCC model

% input: theta --- parameter composed by two symmetric coefficient matrix
% each one has k(k+1)/2 parameters, so theta has k(k+1) parameters
 %       V  ---   matrix of returns
 %       X -----  matrix of realized mesure ,as RK
[k,~,T] = size(X);
numparams = k*k;
B = reshape(theta(1:numparams),k,k);
A = reshape(theta(numparams+(1:numparams)),k,k);
 

P = cov(Y);% Y is lret, P is COVARIANCE MATRIX FOR STANDARDIZED RESIDUALS
M = mean(X,3); % unconditional means of Realized kernel
% k = P\M;  % used in targeting parameter

const = (P-B*P*B'-A*M*A');
H = zeros(k,k,T);
H(:,:,1)=P; % we set the unconditional covariance matrix of daily(close to close) return as starting point of revolved conditional covariance
log_likelihood = zeros(T,1);

for t=2:T
  H(:,:,t) = const + B*H(:,:,t-1)*B'+ A*X(:,:,t-1)*A' ; %conditional covariance MATRIX, insure positive semidefinite
end

% Compute minus log-likelihood
for t=1:T
  log_likelihood(t) = 0.5*log(det(H(:,:,t))) +0.5*trace(H(:,:,t)^(-1)*(Y(t,:)'*Y(t,:)));% Wishart distribution
end

L = sum(log_likelihood);