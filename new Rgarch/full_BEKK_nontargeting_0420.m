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
shape=eye(K);
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
[theta_hat, log_likelihood] = fmincon(@(theta)rbekk_likelihood_nontargeting(theta,lret,logRK), theta, [],[],[],[],LB,UB,@rbekk_constraint_nontargeting,options); % constraint on rbekk
warning('on') %#ok<*WNON>
% extract the coefficient and convert to matrix form
numParams = K*K;
k2 = K*(K+1)/2;
tempP = theta_hat(k2+(1:numParams));
B = reshape(tempP,K,K);SmoothMat = B*B'
k3 = k2+numParams;
tempP = theta_hat(k3+(1:numParams));
A = reshape(tempP,K,K);InovationMat = A*A'


vechC = theta(1:k2);
C=zeros(K);
% Use a logical trick to inverse the vech
ind=tril(true(K));
C(ind)=vechC;
C = C*C'
[L, logH] = rbekk_likelihood_nontargeting(theta_hat,  lret, logRK);
L = -L;
% In estimation  Q(t) = { Shat*(1-alpha-beta) + alpha*y(t-1,:)'*y(t-1,:) + beta*Q{t-1} };
%   R = sqrt(diag(diag(Q{t})));
%   Gamma(t) = { inv(R)*Q{t}*inv(R) };

savefile = 'logH.mat';
save(savefile, 'logH');
savefile = 'logRK.mat';
save(savefile, 'logRK');