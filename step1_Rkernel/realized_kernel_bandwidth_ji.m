function bandwidth = realized_kernel_bandwidth_ji(noiseVariance, IQEstimate, options)


% Revised for Nonflatparzen kernel.

% Get the size of the filtered data 
nMax = options.filteredN;

xiSquared = noiseVariance/sqrt(IQEstimate);
cStar = ((12)^2/0.269)^(1/5);
bandwidth = cStar * (xiSquared)^(2/5) * nMax^(3/5);


end