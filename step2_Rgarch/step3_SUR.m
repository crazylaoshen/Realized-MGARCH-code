clc;
clear;
load('H.mat');
load('RK.mat');
% RK(:,:,all(all(RK==0)==1,2))=[];
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
x = zeros(T,1+K*(K+1)/2);
x(:,1) = ones(T,1);
x(:,2:end) = h_vech';
x = x(2:end,:);
RK_vech = RK_vech(:,2:end);
for s = 1:(K*(K+1)/2);
    

    
    y = RK_vech(s,:)';
    betaHat(:,s) = (x'*x)\x'*y;

end

tauhat = betaHat';
% betaHat = betaHat(:,2:end);