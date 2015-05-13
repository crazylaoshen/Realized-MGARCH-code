% Compute variance of estimates

% clc;clear;close;
% load the data
load('RK.mat');
load('lret.mat');
[K,~,T]= size(RK);
logRK = zeros(K,K,T);
for t=1:T
logRK(:,:,t)=logm(RK(:,:,t));
end
clear RK;
load('A.mat');
load('B.mat');
load('C.mat');
K=10;
C=chol(C)';
ind = tril(true(K));
C = C(ind);


theta_hat = [C;sqrt(B);sqrt(A)];

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
std_errors = sqrt(diag((I\J/I)/T));



disp('Theta_hat  |  std_errors')
[theta_hat, std_errors]
savefile = 'std.mat';
save(savefile, 'std_errors');


   