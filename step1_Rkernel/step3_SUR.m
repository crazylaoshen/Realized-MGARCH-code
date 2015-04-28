clc;
clear;
load('H.mat');
load('RK.mat');
[S,~,T] = size(Sigma);
h_vech = zeros(S*(S+1)/2,T);
RK_vech = zeros(S*(S+1)/2,T);
ind = tril(true(S));
for t = 1:T;
sliceSigma = Sigma(:,:,t);
h_vech(:,t) = sliceSigma(ind);
sliceRK = RK(:,:,t+1);
RK_vech(:,t) = sliceRK(ind);
end

betaHat = zeros(1+S*(S+1)/2,S*(S+1)/2);
x = zeros(T,1+S*(S+1)/2);
x(:,1) = ones(T,1);
x(:,2:end) = h_vech';
x = x(2:end,:);
RK_vech = RK_vech(:,2:end);
for s = 1:(S*(S+1)/2);
    

    
    y = RK_vech(s,:)';
    betaHat(:,s) = (x'*x)\x'*y;

end

% betaHat = betaHat';
% betaHat = betaHat(:,2:end);