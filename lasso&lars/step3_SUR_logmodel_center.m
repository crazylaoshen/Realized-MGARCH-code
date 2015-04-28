clc;
clear;
load('logH.mat');
load('logRK.mat');

% 
% [K,~,T]= size(RK);
% logRK = zeros(K,K,T);
% for t=1:T
% logRK(:,:,t)=logm(RK(:,:,t));
% end

[K,~,T] = size(logH);
h_vech = zeros(K*(K+1)/2,T);
RK_vech = zeros(K*(K+1)/2,T);
ind = tril(true(K));
for t = 1:T;
sliceSigma = logH(:,:,t);
h_vech(:,t) = sliceSigma(ind);
sliceRK = logRK(:,:,t);
RK_vech(:,t) = sliceRK(ind);
end

betaHat = zeros(K*(K+1)/2,K*(K+1)/2);
pvalue = zeros(K*(K+1)/2,K*(K+1)/2);

X = h_vech';

for i = 1:size(X,2)
X(:,i) = X(:,i) - mean(X(:,i));
end




for s = 1:(K*(K+1)/2);
    y = RK_vech(s,:)';    
    y = y-mean(y);  
    betaHat(:,s) = (X'*X)\X'*y;
    uHat = y-X*betaHat(:,s);
    sigma2Hat = (uHat'*uHat )/(T-size(X,2));
    varBetaHat = sigma2Hat*inv(X'*X);
    stdError = sqrt(diag(varBetaHat));
    tStat = betaHat(:,s)./stdError;
    pvalue(:,s) = 1-tcdf(tStat,(T-size(X,2)));
    Rsquare(s) = 1- (uHat'*uHat)/(y'*y);
    

end

tauOLS = betaHat';
pmatrix = pvalue';
% betaHat = betaHat(:,2:end);