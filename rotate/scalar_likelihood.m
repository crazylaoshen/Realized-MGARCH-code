function [ll,lls,logH]=scalar_likelihood(parameters,Y,X)
% Log likelihood for SCALAR_VT_VECH(P,Q) estimation



[K,~,T]=size(X);

beta=parameters(1);
alpha=parameters(2);
%Compute the constant

M = mean(X,3); % unconditional means of Realized kernel
C=cov(Y);% unconditional covariance matrix is just outproduct of returns

intercept=(C-alpha*M-beta*C);
% const=C;
%Initialize the covariance

H=zeros(K,K,T);

%Initialize the log likelihood
lls=zeros(T,1);

%Perform the recursion
H(:,:,1)=C;
for t=2:T;
   
H(:,:,t)=intercept+beta*H(:,:,t-1)+alpha*X(:,:,t-1);
   
end
% Replace these lines to make it work better with poorly conditioned covairance
% likelihoods(i)=likconst+(log(det(Ht(:,:,t)))+data(i,:)*Ht(:,:,t)^(-1)*data(i,:)');
% This is a trick to ensure better numerical stability

for t=1:T
 
    Q=sqrt(diag(H(:,:,t)));
    R=H(:,:,t)./(Q*Q');
    %stdresid=X(:,:,t)./(Q*Q');
    stdresid=Y(t,:)'*Y(t,:)./(Q*Q');
    lls(t)=0.5*(2*sum(log(Q))+log(det(R))+trace(R^(-1)*stdresid));
% lls(t) = 0.5*log(det(H(:,:,t))) +0.5*sum(diag(H(:,:,t)\(Y(t,:)'*Y(t,:))));% Wishart distribution
end
 
ll = sum(lls);