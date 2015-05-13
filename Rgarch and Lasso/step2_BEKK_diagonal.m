clear all
close all;
clc
% load the data
load('RK.mat');
load('lret.mat');
[K,~,T]= size(RK);
logRK = zeros(K,K,T);
for t=1:T
logRK(:,:,t)=logm(RK(:,:,t));
end
clear RK;

% M = mean(RK,3);
P = cov(lret);

% starting value
% startingVals = rvech_starting_values(lret,RK);
 startingVals = [0.6 0.3];

C = P*(1-startingVals(1)-startingVals(2));
C=chol(C)';
ind = tril(true(K));
C = C(ind);
shape=ones(K,1);
theta = [];
for i=1:2
        temp = sqrt(startingVals(i));
        temp = temp * shape;
        theta = [theta;temp(:)]; %#ok<AGROW>
end
theta = [C;theta];
% OPTIONS = optimset('MaxIter', 200,'MaxFunEvals', 1e5, 'Display', 'iter');
options = optimset('MaxIter', 200,'MaxFunEvals', 1e5, 'Display', 'off','Algorithm', 'interior-point');

% Estimate the bekk model
% theta = 0.1*ones(k*k*2,1);
UB = .99998 * ones(size(theta));
LB = -UB;
warning('off') %#ok<*WNOFF>
[theta_hat, log_likelihood] = fmincon(@(theta)rbekk_likelihood_diagonal(theta,lret,logRK), theta, [],[],[],[],LB,UB,@constraint_diagonal,options); % constraint on rbekk
warning('off') %#ok<*WNON>
% extract the coefficient and convert to matrix form
numParams = K;
k2 = K*(K+1)/2;
B = diag(theta_hat(k2+(1:numParams))); % print the estimated parameter
SmoothMat = B*B'
k3 = k2+numParams;
A = diag(theta_hat(k3+(1:numParams)));
InovationMat = A*A'
vechC = theta(1:k2);
C=zeros(K);
% Use a logical trick to inverse the vech
ind=tril(true(K));
C(ind)=vechC;
intercept = C*C'
[L,m , logH] = rbekk_likelihood_diagonal(theta_hat,  lret, logRK);
L = -L;
% In estimation  Q(t) = { Shat*(1-alpha-beta) + alpha*y(t-1,:)'*y(t-1,:) + beta*Q{t-1} };
%   R = sqrt(diag(diag(Q{t})));
%   Gamma(t) = { inv(R)*Q{t}*inv(R) };
savefile = 'theta_hat.mat';
save(savefile, 'logH');
savefile = 'logH.mat';
save(savefile, 'logH');
savefile = 'logRK.mat';
save(savefile, 'logRK');
savefile = 'SmoothMat.mat';
save(savefile, 'SmoothMat');
savefile = 'InovationMat.mat';
save(savefile, 'InovationMat');
savefile = 'intercept.mat';
save(savefile, 'intercept');
%% calculate the SE
p = length(theta_hat);
I = -hessian(@(theta)rbekk_likelihood_diagonal(theta,lret,logRK),theta_hat)/T;
savefile = 'I.mat';
save(savefile, 'I');
first_diff = zeros(T,p);
parfor i=1:p
  epsilon = zeros(p,1);
  epsilon(i) = 1e-7;
  [~,tutu1] = rbekk_likelihood_diagonal(theta_hat+epsilon,lret,logRK);
  [~,tutu2] = rbekk_likelihood_diagonal(theta_hat-epsilon,lret,logRK);
  first_diff(:,i) = (tutu1-tutu2) / (2*epsilon(i));
end
J = (first_diff'*first_diff)/T;
savefile = 'J.mat';
save(savefile, 'J');
se = sqrt(diag((inv(I)*J*inv(I))/T));



disp('Theta_hat  |  std_errors')
[theta_hat, se]
savefile = 'se.mat';
save(savefile, 'se');
vechC_se = se(1:k2);
C_se=zeros(K);
% Use a logical trick to inverse the vech
ind=tril(true(K));
C_se(ind)=vechC_se;
tempP = se(k2+(1:numParams));
B_se = diag(tempP);
k2 = k2+numParams;
tempP = se(k2+(1:numParams));
A_se = diag(tempP);
savefile = 'B_se.mat';
save(savefile, 'B_se');
savefile = 'A_se.mat';
save(savefile, 'A_se');
savefile = 'C_se';
save(savefile, 'C_se');