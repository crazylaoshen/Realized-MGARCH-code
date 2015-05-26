clc;
clear;
load('H.mat');
load('logRK.mat');

[K,~,T]= size(H);
logH = zeros(K,K,T);
for t=1:T
logH(:,:,t)=logm(H(:,:,t));
end
clear H;

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
AICse = zeros(K*(K+1)/2,K*(K+1)/2);

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
% betastar = lassoBICse(G,y);
    yfit = G*betastar;
    resid = y - yfit;
    
    AICse(:,s) = std(bootstrp(...
         100,@(bootr)lassoBICse(G,yfit+bootr),resid))'.*(betaOLS(:,s)).^2;
    
    
    betaHat(:,s) = betastar.*abs(betaOLS(:,s));
end
toc
tauLARS = betaHat';
tauLARS = tauLARS./(ones(p,1)*d);
tauSE = AICse';
tauSE = tauSE./(ones(p,1)*d.^2);
tauOLS = betaOLS'./(ones(p,1)*d);




% Overall Percent of element being non-zero
pctAll = sum(sum(tauLARS~=0))/numel(tauLARS);

%Row1: Percent of element being non-zero
row1 = tauLARS(1,:)';
tuRow1 = sum(row1~=0);
pctRow1 = sum(tuRow1)/length(row1);

savefile = 'AICse.mat';
save(savefile, 'tauSE');


'Diagonal element of tauhat '
tu = diag(tauLARS)'
'Mean of diagonal element'
mean(tu)
% Diagonal :Percent of diagonal element being non-zero'
pctDiag = sum(tu~=0)/length(tu);

'Percent of element being non-zero'
'||Overall ||  Row1  || Diag||'
[pctAll, pctRow1, pctDiag]     

% SE of row 1
se1 = tauSE(1,:)';
% SE of diagonal
sed = diag(tauSE);





