function [L, h] = garch(theta, y, omega)
% code from Dr. Denis Pelletier's course in asset pricing @ NCSU

% Minus the log-likelihood of a GARCH(1,1) model

alpha = theta(1);
beta = theta(2);

T = size(y,1);
log_likelihood = zeros(T,1);
u = zeros(T,1);
h = zeros(T,1);

h(1) = omega;
u(1) = y(1) / sqrt(h(1));
log_likelihood(1) = 0.5 * log(2*pi) + 0.5 * log(h(1)) + 0.5 * u(1)^2;

for t=2:T
  h(t) = omega*(1-alpha-beta) + alpha * y(t-1)^2 + beta * h(t-1);
  u(t) = y(t) / sqrt(h(t));
  log_likelihood(t) = 0.5 * log(2*pi) + 0.5 * log(h(t)) + 0.5 * u(t)^2;
end

L = sum(log_likelihood);
