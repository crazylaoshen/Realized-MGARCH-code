

%This code performs a grid search of potential coefficient beta and alpha
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
a=linspace(.01,.99,50)';
b=linspace(.01,.99,50)';
[A,B] = meshgrid(a,b);
c=cat(2,A',B');
params=reshape(c,[],2);

params = unique(params,'rows');

ll=zeros(length(params),1);
%Evaluate the likelihood at each parameter
parfor i=1:size(params,1)
    parameters=params(i,:);
    ll(i)=scalar_likelihood_logmodel(parameters,lret,logRK);
end

ll(imag(ll)~=0)=nan;
[sortedll,pos]=sort(ll);
output_parameters=params(pos,:);
bestparams=output_parameters(1,:); 
 % beta = 0.51  alpha = 0.3 for 4 stocks
 % 0.48  0.26 for 2 stocks
 % 0.7   0.16 for realized covariance of 4 stocks
 % 0.53  0.29 for realized covariance of 5 min and 1min subsampling
 % 0.54 0.24 for 3 banks
 %0.57 0.3 for new 4 stocks
 %0.47 0.31 for new 3 banks
 %0.98 0.49 for simulation
[L ,~ , logH] = scalar_likelihood_logmodel(bestparams,lret,logRK);
L = -L;
savefile = 'logH.mat';
save(savefile, 'logH');
savefile = 'logRK.mat';
save(savefile, 'logRK');