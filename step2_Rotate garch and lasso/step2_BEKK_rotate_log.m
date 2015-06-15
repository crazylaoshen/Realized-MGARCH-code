clear all
close all;
clc
% pctRunOnAll warning('off','all') 
% load the data
load('logRK.mat');
load('lret.mat');
[K,~,T]= size(logRK);
logH = zeros(K,K,T);
RRK = zeros(K,K,T);
M = mean(logRK,3);
MM = spectral(M); % M^(1/2)  spectral decompostion of M

% Rotated Realized kernel
for t=1:T
RRK(:,:,t)=MM\logRK(:,:,t)/MM';
end
P = cov(lret);

PP = spectral(P); % P^(1/2)  spectral decompostion of P
e = lret/PP;


% starting value
 startingVals = [0.6 0.3];


shape=ones(K,1);
theta = [];
for i=1:2
        temp = sqrt(startingVals(i));
        temp = temp * shape;
        theta = [theta;temp(:)]; %#ok<AGROW>
end

% OPTIONS = optimset('MaxIter', 200,'MaxFunEvals', 1e5, 'Display', 'iter');
options = optimset('MaxIter', 200,'MaxFunEvals', 1e5, 'Display', 'off','Algorithm', 'interior-point');

% Estimate the bekk model
% theta = 0.1*ones(k*k*2,1);
UB = .99998 * ones(size(theta));
LB = -UB;
warning('off') %#ok<*WNOFF>
tic
[theta_hat, log_likelihood] = fmincon(@(theta)rbekk_likelihood_rotate_log(theta,e,RRK), theta, [],[],[],[],LB,UB,@constraint_rotate,options); % constraint on rbekk
warning('off') %#ok<*WNON>
toc
% extract the coefficient and convert to matrix form
numParams = K;

B = diag(theta_hat(1:numParams)); % print the estimated parameter
B_true = PP*B*inv(PP);
SmoothMat = B*B'
A = diag(theta_hat(K+(1:numParams)));
A_true = PP*A*inv(MM);
InovationMat = A*A'
intercept = eye(K)-InovationMat-SmoothMat;

[L, m, logRH] = rbekk_likelihood_rotate(theta_hat,e,RRK);
L = -L;

for t=1:T
logH(:,:,t)=PP*logRH(:,:,t)*PP';
end

savefile = 'theta_hat.mat';
save(savefile, 'theta_hat');
savefile = 'H.mat';
save(savefile, 'H');
savefile = 'B_tilde.mat';
save(savefile, 'B');
savefile = 'B.mat';
save(savefile, 'B_true');
savefile = 'A_tilde.mat';
save(savefile, 'A');
savefile = 'A.mat';
save(savefile, 'A_true');
savefile = 'intercept.mat';
save(savefile, 'intercept');
%% calculate the SE
p = length(theta_hat);
I = -hessian(@(theta)rbekk_likelihood_rotate(theta,e,RRK),theta_hat)/T;
savefile = 'I.mat';
save(savefile, 'I');
first_diff = zeros(T,p);
parfor i=1:p
  epsilon = zeros(p,1);
  epsilon(i) = 1e-7;
  [~,tutu1] = rbekk_likelihood_rotate(theta_hat+epsilon,e,RRK);
  [~,tutu2] = rbekk_likelihood_rotate(theta_hat-epsilon,e,RRK);
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
tempP = se(1:numParams);
B_se = diag(tempP);
tempP = se(K+(1:numParams));
A_se = diag(tempP);
savefile = 'B_se.mat';
save(savefile, 'B_se');
savefile = 'A_se.mat';
save(savefile, 'A_se');
