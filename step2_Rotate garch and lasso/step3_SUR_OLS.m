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
se = zeros(K*(K+1)/2,K*(K+1)/2);
X = h_vech';

tic
for s = 1:(K*(K+1)/2);
    y = RK_vech(s,:)';
    y = y-mean(y);
    betaOLS(:,s) = (X'*X)\X'*y;
    residOLS = y - X*betaOLS(:,s);
    sigmaHat = (residOLS'*residOLS)/(T-K);
    se(:,s) = sqrt(diag(sigmaHat*inv(X'*X)));
    
    
end
toc
tauLARS = betaHat';
tause = se';
se1 = se(1,:)';
sed = diag(tause);





