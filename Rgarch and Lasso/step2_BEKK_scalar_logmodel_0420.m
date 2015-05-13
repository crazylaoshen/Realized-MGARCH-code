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

 theta = [0.6 0.3];

% theta = [C;theta];
% OPTIONS = optimset('MaxIter', 200,'MaxFunEvals', 1e5, 'Display', 'iter');
options = optimset('MaxIter', 200,'MaxFunEvals', 1e5, 'Display', 'off','Algorithm', 'interior-point');

% Estimate the bekk model
% theta = 0.1*ones(k*k*2,1);
UB = .99998 * ones(size(theta));
LB = -UB;
warning('off') %#ok<*WNOFF>
[theta_hat, log_likelihood] = fmincon(@(theta)scalar_likelihood_logmodel(theta,lret,logRK), theta, [],[],[],[],LB,UB,[],options); % constraint on rbekk
warning('on') %#ok<*WNON>
% extract the coefficient and convert to matrix form
B = theta_hat(1); % print the estimated parameter
SmoothMat = B
A = theta_hat(2);
InovationMat = A
% Use a logical trick to inverse the vech

[L,~, logH] = scalar_likelihood_logmodel(theta_hat,  lret, logRK);
L = -L;
% In estimation  Q(t) = { Shat*(1-alpha-beta) + alpha*y(t-1,:)'*y(t-1,:) + beta*Q{t-1} };
%   R = sqrt(diag(diag(Q{t})));
%   Gamma(t) = { inv(R)*Q{t}*inv(R) };

savefile = 'logH.mat';
save(savefile, 'logH');