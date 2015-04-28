clear all
close all;
clc
% load the data
% load('cleaned_data.mat');
load('RK.mat');
load('lret.mat');
% RK(:,:,all(all(RK==0)==1,2))=[];
% closePrice(all(closePrice==0,2),:)=[];


% lret = diff(log(closePrice));
% lret = lret-ones(size(lret,1),1)*mean(lret);
% test = mean(lret);

[T,k] = size(lret);
M = mean(RK,3);
P = cov(lret);
% Rotated RK for target covariance constraint. see page9 of Multivariate
% HEAVY paper
% for t=1:T
% RK(:,:,t) = chol(M)'^(-1)*chol(P)'*RK(:,:,t)*chol(P)*chol(M)^(-1);
% end
% starting value
startingVals = rvech_starting_values(lret,RK);
% startingVals = [0.2 0.6];
shape=eye(k);
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
[theta_hat, log_likelihood] = fmincon(@(theta)rbekk_likelihood_targeting(theta,lret,RK), theta, [],[],[],[],LB,UB,@rbekk_constraint_targeting,options); % constraint on rbekk
warning('on') %#ok<*WNON>
numparams = k*k;
B = reshape(theta_hat(1:numparams),k,k);
SmoothMat = B*B'

A = reshape(theta_hat(numparams+(1:numparams)),k,k); % print the estimated parameter
% A = A*chol(P)'^(-1)*chol(M)'^(-1);
InovationMat = A*A'
C = (P-B*P*B'-A*M*A')
[L, H] = rbekk_likelihood_targeting(theta_hat,  lret, RK);
L = -L;
% In estimation  Q(t) = { Shat*(1-alpha-beta) + alpha*y(t-1,:)'*y(t-1,:) + beta*Q{t-1} };
%   R = sqrt(diag(diag(Q{t})));
%   Gamma(t) = { inv(R)*Q{t}*inv(R) };

savefile = 'H.mat';
save(savefile, 'H');