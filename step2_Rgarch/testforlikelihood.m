clc;clear;
load('RK.mat');
load('lret.mat');
X = RK;
Y = lret;
[k,~,T] = size(X);
numparams = k*k;
startingVals = rvech_starting_values(lret,RK);
% startingVals = [0.6 0.2];
shape=eye(k);
theta = [];
for i=1:2
        temp = sqrt(startingVals(i));
        temp = temp * shape;
        theta = [theta;temp(:)]; %#ok<AGROW>
end
B = reshape(theta(1:numparams),k,k);
A = reshape(theta(numparams+(1:numparams)),k,k);
 

P = cov(Y);% Y is lret, P is COVARIANCE MATRIX FOR STANDARDIZED RESIDUALS
M = mean(X,3); % unconditional means of Realized kernel
% k = P\M;  % used in targeting parameter

const = (P-B*P*B'-A*M*A');
H = zeros(k,k,T);
H(:,:,1)=P; % we set the unconditional covariance matrix of daily(open to close) return as starting point of revolved conditional covariance
log_likelihood = zeros(T,1);

for t=2:T
  H(:,:,t) = const + B*H(:,:,t-1)*B'+ A*X(:,:,t-1)*A' ; %conditional covariance MATRIX, insure positive semidefinite
end

% Compute minus log-likelihood
for t=1:T
  log_likelihood(t) = 0.5*log(100*det(H(:,:,t))) +0.5*trace(H(:,:,t)^(-1)*(Y(t,:)'*Y(t,:)));% Wishart distribution
end
L = sum(log_likelihood);
