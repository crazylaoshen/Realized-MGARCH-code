clc
clear
% 
% H = zeros(T,K)
% for t = 1:T
%     
% H(1,:) = diag(H(:,:,T))'
% end;
load('RK.mat');
load('lret.mat');
[K,~,T]= size(RK);
logRK = zeros(K,K,T);
for t=1:T
logRK(:,:,t)=logm(RK(:,:,t));
end

% for t=1:T
%     H(:,:,t) = expm(logH(:,:,t));
% end
V = zeros(T,K);
for t = 1:T
    
V(t,:) = sqrt(diag(logRK(:,:,t))');
end;

% hist(V(:,1),100)
skewness(V)
kurtosis(V)