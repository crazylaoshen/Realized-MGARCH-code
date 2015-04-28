function [L, h] = rgarch(theta, y, x)
% Minus the log-likelihood of the log likelihood model
% y = return ( or residual of VAR)
% x = realized variances


omega   = theta(1);
beta    = theta(2);
alpha   = theta(3);

% omegaH = mean(x)*(1-alpha-beta);

T = size(y,1);
h = zeros(T,1);

h(1) =omega/(1-beta-alpha); %initial value of h_1

for t=2:T
  h(t) = omega + beta * (h(t-1))+  alpha * x(t-1) ;
end
u = y./ sqrt(h);
log_likelihood =  0.5*(log(2*pi) + log(h) + u.^2 );
L = sum(log_likelihood);



