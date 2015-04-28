
function  jitterLags = realized_kernel_jitter_lag_length_ji(noiseEstimate,iqEstimate,N)

% Computes the optimal amount of end point jitter given   nonflatpazen kernel, noise estimate, integrated
% quarticity estimate and number of observations.





% n^1/5 non flat top
cStar = ((12)^2/0.269)^(1/5);  % for nonflatparzen kernel
k00 = .269;
avar = 5 * cStar * k00 * noiseEstimate^(4/5) * iqEstimate^(4/5);
power = 1/5;
    
    
% Problem is to min 8*noise^2*jitterLag^(-2) + avar * (N-jitterLag)^(-2*power)
% jittering over m data points, N observations in total
    
    MSE = inf*ones(N-1,1);
for m=1:(N-1)
    MSE(m) = 8 * noiseEstimate^2 * m^(-2) + avar * (N-m)^(-2*power);
end

% Find the minimum, which may be 1
    [temp,jitterLags] = min(MSE);
    
    
    
    
    
    
end