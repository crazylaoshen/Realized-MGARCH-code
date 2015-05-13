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
betaOLS = zeros(K*(K+1)/2,K*(K+1)/2);

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
    betaOLS(:,s) = (X'*X)\X'*y;
    G = X.*(ones(T,1)*abs(betaOLS(:,s)'));
    [beta, info] = lassolar(G, y);
    [bestAIC, bestIdx] = min(info.AIC);
%     ind=beta(:,bestIdx)~=0;
%     aa = sum(ind);
%     XA = X(:,ind);
    betastar = beta(:,bestIdx);
    betaHat(:,s) = betastar.*abs(betaOLS(:,s));
    

end
toc
tauLARS = betaHat';
tauLARS = tauLARS./(ones(p,1)*d);

tauOLS = betaOLS'./(ones(p,1)*d);
% pmatrix = pvalue';
tutu=sum(tauLARS~=0);
'Percent of element being non-zero'
pctzeros = sum(tutu)/numel(tauLARS)
% savefile = 'tauLARS.mat';
% save(savefile, 'tauLARS');
% savefile = 'pvalue.mat';
% save(savefile, 'pmatrix');

'Diagonal element of tauhat '
tu = diag(tauLARS)'
'Mean of diagonal element'
mean(tu)
'Percent of diagonal element being non-zero'
sum(tu~=0)/length(tu)






