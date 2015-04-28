function [b, info] = lassolar(X, y)



%% LARS variable setup
[n p] = size(X);

  Gram = X'*X;


 %% LARS algorithm setup
 maxSteps = 8*p;% Maximum number of algorithm steps
 b = zeros(p, 2*p);% set up the LASSO coefficient vector

mu = zeros(n, 1); % current "position" as LARS travels towards lsq solution
 I = 1:p; % inactive set
A = []; % active set
 lassoCond = 0; % LASSO condition boolean
 step = 1; % step count
%% LARS main loop
while length(A) < p &&  step < maxSteps
  r = y - mu; 
 
  % find max correlation
  c = X'*r;
  [cmax, cidxI] = max(abs(c(I)));
  cidx = I(cidxI); % index of next active variable
 if ~lassoCond 
     A = [A cidx]; % add to active set
    I(cidxI) = []; % ...and drop from inactive set
  else
    % if a variable has been dropped, do one step with this
    % configuration (don't add new one right away) 
    lassoCond = 0;
 end   
    b_OLS = Gram(A,A)\(X(:,A)'*y);  % partial OLS solution
    
    %direction from current position to the OLS solution of X_A
    d = X(:,A)*b_OLS - mu;
  
  % compute length of walk along equiangular direction
  gamma_tilde = b(A(1:end-1),step)./(b(A(1:end-1),step) - b_OLS(1:end-1,1));
 gamma_tilde(gamma_tilde <= 0) = inf;
  [gamma_tilde, dropIdx] = min(gamma_tilde);
 
  if isempty(I)
    % if all variables active, go all the way to the OLS solution
    gamma = 1;
  else
    cd = X'*d;
    temp = [ (c(I) - cmax)./(cd(I) - cmax); (c(I) + cmax)./(cd(I) + cmax) ];
    temp = sort(temp(temp > 0)); % faster than min(temp(temp > 0)) (!)
     gamma = temp(1);
  end
    
    % check if variable should be dropped
  if gamma_tilde < gamma,
    lassoCond = 1; % drop it
    gamma = gamma_tilde;
  end
  % check if beta must grow
     if size(b,2) < step + 1
      b = [b zeros(p, size(b,2))];
    end
    b(A,step + 1) = b(A,step) + gamma*(b_OLS - b(A,step)); % update beta
     
  % update position
  mu = mu + gamma*d; % length times direction
  
  % increment step counter
  step = step + 1;
  
  
  if lassoCond
  I = [I A(dropIdx)]; % add dropped variable to inactive set
    A(dropIdx) = []; % ...and remove from active set
  end
end
 
 % trim beta
if size(b,2) > step
  b(:,step + 1:end) = [];
end
% return number of iterations
steps = step - 1;
 

%% Compute auxilliary measures
if nargout == 2 % only compute if asked for
  info.steps = steps;
  b0 = pinv(X)*y; % regression coefficients of low-bias model
  penalty0 = sum(abs(b0)); % L1 constraint size of low-bias
  indices = (1:p)';
  
 % for entire path
    q = info.steps + 1;
    sigma2e = sum((y - X*b0).^2)/n;
    info.lambda = zeros(1,q);
    info.df = zeros(1,q);
    info.Cp = zeros(1,q);
    info.AIC = zeros(1,q);
    info.BIC = zeros(1,q);
    info.s = zeros(1,q);
    for step = 1:q
      A = indices(b(:,step) ~= 0); % active set
      % compute godness of fit measurements Cp, AIC and BIC
      r = y - X(:,A)*b(A,step); % residuals
      rss = sum(r.^2); % residual sum-of-squares
      info.df(step) = length(A);
      info.Cp(step) = rss/sigma2e - n + 2*info.df(step);
      info.AIC(step) = rss + 2*sigma2e*info.df(step);
      info.BIC(step) = rss + log(n)*sigma2e*info.df(step);
      % compute L1 penalty constraints s and lambda
      info.s(step) = sum(abs(b(A,step)))/penalty0;
      if (step == 1)
        info.lambda(step) = max(2*abs(X'*y));
      else
        info.lambda(step) = median(2*abs(X(:,A)'*r));
      end
    end
    
  else % for single solution
    % compute L1 penalty constraints s and lambda at solution
    A = indices(b ~= 0); % active set
    info.s = sum(abs(b))/penalty0;
    info.df = length(A);
    info.steps = steps;
    if isempty(A)
      info.lambda = max(2*abs(X'*y));
    else
      info.lambda = median(2*abs(X(:,A)'*(y - X(:,A)*b(A))));
    end

  
end
