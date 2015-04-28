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
p = K*(K+1)/2;
h_vech = zeros(K*(K+1)/2,T);
RK_vech = zeros(K*(K+1)/2,T);
ind = tril(true(K));
for t = 1:T;
    sliceSigma = logH(:,:,t);
    h_vech(:,t) = sliceSigma(ind);
    sliceRK = logRK(:,:,t);
    RK_vech(:,t) = sliceRK(ind);
    clear sliceRK slickSigma;
 
end

betaHat = zeros(K*(K+1)/2,K*(K+1)/2);
pvalue = zeros(K*(K+1)/2,K*(K+1)/2);

X = h_vech';
% normalize X
X = X - ones(T,1)*mean(X);
d = sqrt(sum(X.^2));
d(d == 0) = 1;
X = X./(ones(T,1)*d);


% lambda_min = 0;
% lambda_max = 6;
% nlambda = 60;
% lambda_vec = linspace(lambda_min, lambda_max, nlambda);
% kfold = 10;
% lambda = zeros(K*(K+1)/2,1);
tic
for s = 1:(K*(K+1)/2);
    y = RK_vech(s,:)';
    y = y-mean(y);

    [beta, info] = lassolar(X, y);
    [bestCp, bestIdx] = min(info.Cp);


    betaHat(:,s) = beta(:,bestIdx);

    
%     uHat = y-X*betaHat(:,s);
%     sigma2Hat = (uHat'*uHat )/(T-size(X,2));
%     varBetaHat = sigma2Hat*inv(X'*X);
%     stdError = sqrt(diag(varBetaHat));
%     tStat = betaHat(:,s)./stdError;
%     pvalue(:,s) = 1-tcdf(tStat,(T-size(X,2)));
%     Rsquare(s) = 1- (uHat'*uHat)/(y'*y);
end
toc
taulars = betaHat';
taulars = taulars./(ones(p,1)*d);
% pmatrix = pvalue';
tutu=sum(taulars==0);
pctzeros = sum(tutu)/numel(taulars)
savefile = 'taulars_nodapt.mat';
save(savefile, 'taulars');
% savefile = 'pvalue.mat';
% save(savefile, 'pmatrix');

'Diagonal element of tauhat '
tu = diag(taulars)'
'Mean of diagonal element'
mean(tu)
'Percent of diagonal element being zero'
sum(tu==0)/size(taulars,1)
% 
% ptu = diag(pmatrix);
% 'Mean of P value of element which donot equal to zero'
% mean(ptu(tu~= 0))
