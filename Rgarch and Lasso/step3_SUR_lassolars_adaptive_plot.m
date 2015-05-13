clc;
clear; close all;
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


s = 1;  % choose the eqation 
    y = RK_vech(s,:)';
    y = y-mean(y);
    betaOLS = (X'*X)\X'*y;
    G = X.*(ones(T,1)*abs(betaOLS'));
    [beta, info] = lassolar(G, y);
    [bestAIC, bestAIdx] = min(info.AIC);
    [bestBIC, bestBIdx] = min(info.BIC);

%     betastar = beta(:,bestAIdx);
%     betaHat(:,s) = betastar.*abs(betaOLS);

beta = beta.*(abs(betaOLS)*ones(1,56));
beta = beta./(d'*ones(1,56));
betaOLS = betaOLS./d';
betaAIC = beta(:,bestAIdx);
betaBIC = beta(:,bestBIdx);
LA = info.s(bestAIdx);
LB = info.s(bestBIdx);

figure
plot(info.s,beta');
hold on
plot([LA LA],[-1 1.5],':m','LineWidth',5);
plot([LB LB],[-1 1.5],'-b','LineWidth',5);
hold off
set(gca,'XDir','rev')
xlabel('Relative bound S')
ylabel('Coefficient of one mesurement equation')
title('Lasso shrinkage of coefficients')      

figure(2)
stairs(info.s,info.df,'-*b');
hold on
plot([LA LA],[0 60],':m','LineWidth',5);
plot([LB LB],[0 60],'-b','LineWidth',5);
hold off
set(gca,'XDir','rev')
xlabel('Relative bound S')
ylabel('Df of coefficients')
title('Degree of freedom shrinkage') 

figure(3)
plot(info.s, info.mse)
hold on
plot([LA LA],[0.2 0.7],':m','LineWidth',5);
plot([LB LB],[0.2 0.7],'-b','LineWidth',5);
hold off
xlabel('Relative bound S')
ylabel('MSE of ALASSO')
title('Mean Squared Error of ALASSO') 

['RSS of OLS  ',  '  RSS of AIC based','  RSS of BIC based']
[info.rss(p)    info.rss(bestAIdx)        info.rss(bestBIdx)]
