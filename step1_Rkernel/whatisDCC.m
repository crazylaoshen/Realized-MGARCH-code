beta = 0.1;
alpha  = 0.1;
X= RK(:,:,2:end);
[N,~,T] = size(X);

log_likelihood = zeros(T,1);
Q = zeros(N,N,T);
Gamma = zeros(N,N,T);
%Q = cell(T,1);  % correlation matrix
%Gamma = cell(T,1); % intercept or long run correlation matrix

% set up initial value
Q(:,:,1) = Shat;% unconditional CORRELATION MATRIX FOR STANDARDIZED RESTURN
% R = sqrt(diag(diag(Q(:,:,1))));
% Gamma(:,:,1) = inv(R)*Q(:,:,1)*inv(R) ;
for t=2:T
  Q(:,:,t) = Shat*(1-alpha-beta) + alpha*X(:,:,t-1) + beta*Q(:,:,t-1);
  R = sqrt(diag(diag(Q(:,:,t))));
  Gamma(:,:,t) =  inv(R)*Q(:,:,t)*inv(R) ;
end

% Compute minus log-likelihood
for t=1:T
  log_likelihood(t) = 0.5*log(det(Q(:,:,t))) +0.5*(V(t,:)*inv(Q(:,:,t))*V(t,:)');% y(t,:)*inv(Gamma(:,:,t))*y(t,:)'
end

L = sum(log_likelihood);