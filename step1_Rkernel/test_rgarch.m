

omega   = 0.3;
beta    = 0.3;
alpha   = 0.3;

% omegaH = mean(x)*(1-alpha-beta);
i=1;
y = a_hat(:,i);
T = size(y,1);
x=rv(2:end,i);
log_likelihood = zeros(T,1);
u = zeros(T,1);
h = zeros(T,1);

h(1) =omega/(1-beta-alpha); %initial value of h_1

for t=2:T
  h(t) = omega + beta * (h(t-1))+ alpha * x(t-1);
end
u(t) = y(t)./ sqrt(h(t));
log_likelihood =  0.5*(log(2*pi) + log(h) + u.^2 );
L = sum(log_likelihood);



