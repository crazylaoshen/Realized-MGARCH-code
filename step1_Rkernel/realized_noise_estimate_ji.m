function [noiseVariance, debiasedNoiseVariance, IQEstimate, noiseEstimateOomen] = realized_noise_estimate_ji(price, time, timeType, options)

noiseVarianceSamplingType = options.noiseVarianceSamplingType;
noiseVarianceSamplingInterval = options.noiseVarianceSamplingInterval;

% Compute the noiseVariance using the debiased version of Bandi-Russell, using all prices
% noiseVariance = realized_variance(price, time, timeType , noiseVarianceSamplingType , noiseVarianceSamplingInterval);


logPrice =log(price);
% Filter prices and compute the RV
filteredLogPrice = realized_price_filter_ji(logPrice,time,timeType,noiseVarianceSamplingType,noiseVarianceSamplingInterval);
returns = diff(filteredLogPrice);
noiseVariance=returns' * returns;


% Compute the effective 'n' to use in the Bandi-Russel estimator
noiseFilteredPrice = realized_price_filter_ji(price, time, timeType , noiseVarianceSamplingType , noiseVarianceSamplingInterval);
n = sum(diff(noiseFilteredPrice)~=0);
noiseVariance = noiseVariance/(2*n);


IQEstimationSamplingType = options.IQEstimationSamplingType;
IQEstimationSamplingInterval = options.IQEstimationSamplingInterval;
% To compute the "optimal" bandwidth, need an estimate of the noise and the QV using low frequency data
% lowFrequencyRealizedVariance = realized_variance(price,time, timeType , IQEstimationSamplingType, IQEstimationSamplingInterval);
filteredLogPrice = realized_price_filter_ji(logPrice,time,timeType,IQEstimationSamplingType,IQEstimationSamplingInterval);
returns = diff(filteredLogPrice);
lowFrequencyRealizedVariance=returns' * returns;



IQEstimate = lowFrequencyRealizedVariance^2;
debiasedNoiseVariance=[];

end