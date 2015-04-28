
function par_autocor = pacf(x,h)


par_autocor = zeros(h,1);
for i=1:h
  lhs = x(1+i:end);
  rhs = ones(length(x)-i,i+1);
  for j=1:i
    rhs(:,1+j) = x(1+i-j:end-j);
  end
  
  beta_hat = inv(rhs'*rhs)*rhs'*lhs;
  par_autocor(i) = beta_hat(1+i);
end
end