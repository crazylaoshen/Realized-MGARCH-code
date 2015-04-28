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
p = K*(K+1)/2; % tauHat is p by p matrix
h_vech = zeros(p,T);
RK_vech = zeros(p,T);
ind = tril(true(K));
for t = 1:T;
    sliceSigma = logH(:,:,t);
    h_vech(:,t) = sliceSigma(ind);
    sliceRK = logRK(:,:,t);
    RK_vech(:,t) = sliceRK(ind);
    clear sliceRK slickSigma;
 
end

betaHat = zeros(p,p);
pvalue = zeros(p,p);

X = h_vech';
X = X - ones(T,1)*mean(X);

lambda_min = 0;
lambda_max = 5;
nlambda = 51;
lambda_vec = linspace(lambda_min, lambda_max, nlambda);
kfold = 10;
% lambda = zeros(p,1);


parfor s = 1:p;
    y = RK_vech(s,:)';
    y = y-mean(y);
    betaOLS = (X'*X)\X'*y;
    G = X.*(ones(T,1)*abs(betaOLS'));
    lambda = LassoLambda( y, G, kfold, lambda_vec );
    betastar = Lassoshoot( y, G, lambda );
    betaHat(:,s) = betastar.*abs(betaOLS);
%     uHat = Y-X*betaHat(:,s);
%     sigma2Hat = (uHat'*uHat )/(T-size(X,2));
%     varBetaHat = sigma2Hat*inv(X'*X);
%     stdError = sqrt(diag(varBetaHat));
%     tStat = betaHat(:,s)./stdError;
%     pvalue(:,s) = 1-tcdf(tStat,(T-size(X,2)));
%     Rsquare(s) = 1- (uHat'*uHat)/(Y'*Y);
s
end

tauHat = betaHat';

% pmatrix = pvalue';
tutu=sum(tauHat==0);
pctzeros = sum(tutu)/p^2
savefile = 'tauHat.mat';
save(savefile, 'tauHat');
% savefile = 'pvalue.mat';
% save(savefile, 'pmatrix');

'Diagonal element of tauhat '
tu = diag(tauHat)'
'Mean of diagonal element'
mean(tu)
'Percent of diagonal element being zero'
sum(tu==0)/p

% ptu = diag(pmatrix);
% ptu(tu==0)
% 'Mean of P value of element which donot equal to zero'
% mean(ptu(tu~= 0))
