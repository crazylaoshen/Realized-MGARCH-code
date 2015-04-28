clc;
clear;
load('H.mat');
load('RK.mat');
[K,~,T] = size(H);
h_vech = zeros(K*(K+1)/2,T);
RK_vech = zeros(K*(K+1)/2,T);
ind = tril(true(K));
for t = 1:T;
sliceSigma = H(:,:,t);
h_vech(:,t) = sliceSigma(ind);
sliceRK = RK(:,:,t);
RK_vech(:,t) = sliceRK(ind);
end

betaHat = zeros(1+K*(K+1)/2,K*(K+1)/2);
pvalue = zeros(1+K*(K+1)/2,K*(K+1)/2);

x = h_vech';

whichstats = {'tstat','rsquare'};



for s = 1:(K*(K+1)/2);
    y = RK_vech(s,:)';    
    stats = regstats(y,x,'linear',whichstats);
     Rsquare(s) = stats.rsquare;

    betaHat(:,s) = stats.tstat.beta;
    pvalue(:,s) = stats.tstat.pval;
    

end

tauhat = betaHat';
pmatrix = pvalue';
% betaHat = betaHat(:,2:end);