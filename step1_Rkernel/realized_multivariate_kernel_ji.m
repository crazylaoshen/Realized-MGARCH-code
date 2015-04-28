function [rmk,diagnostics]=realized_multivariate_kernel_ji(varargin)
% Computes the multivariate version of the realized kernel of BNHLS using the non-flat-top kernels
% and end-point jittering
%
% USAGE:
%   [RMK,DIAGNOSTICS] = realized_multivariate_kernel(PRICE01,TIME01,PRICE02,TIME02,...,PRICEN,TIMEN,..
%                                                    TIMETYPE,SAMPLINGINTERVAL,OPTIONS)
%   [RMK,DIAGNOSTICS] = realized_multivariate_kernel(PRICES,TIMES,TIMETYPE,SAMPLINGINTERVAL,OPTIONS)
% INPUTS:
%   PRICE01          - m01 by 1 vector of high frequency prices
%   TIME01           - m01 by 1 vector of times where TIME01(i) corresponds to PRICE01(i)
%   PRICE02          - m02 by 1 vector of high frequency prices
%   TIME02           - m02 by 1 vector of times where TIME02(i) corresponds to PRICE02(i)
%   [N price-time pairing as required, where N can be any integer]
%   PRICEN           - mN by 1 vector of high frequency prices
%   TIMEN            - mN by 1 vector of times where TIMEN(i) corresponds to PRICEN(i)
%     -or-
%   PRICES           - Cell array containing prices, e.g. PRICES{1}=PRICE01,
%                        PRICES{2}=PRICE02, etc.
%   TIMES            - Cell array containing observation times, e.g. TIMES{1}=TIME01,
%                        TIMES{2}=TIME02, etc.
%   TIMETYPE         - String describing the way times are measured
%                       'wall'    24-hour clock of the form HHMMSS, e.g. 101543 or 153217
%                       'seconds' Time measured in seconds past midnight on the first day.
%                       'unit'  Unit normalized date format, e.g. .1, .234, .9
%                         Unit normalized times are more general than the other types and can be
%                         applied to data from more than one calendar day
%   SAMPLINGINTERVAL - Scalar integer or n by 1 vector which indicates the sampling interval of the
%                        refresh-time synchronized prices.  This value is the multivariate extension
%                        of the SAMPLINGINTERVAL for 'BusinessTime' sampling in the univariate kernel.
%                        1 indicates to use every price.  2 indicates to use every other price, and
%                        so on. 1 is usually the correct choice.
%   OPTIONS          - [OPTIONAL] Realized Kernel option structure initialized by calling
%                        realized_options('Kernel'). See help realized_options for a description of
%                        available options.
%
% OUTPUTS:
%   RMK              - Realized covariance estimate
%   DIAGNOSTICS      - Structure of useful diagnostic information.  Fields are
%                     KERNEL                - The kernel used in the RMK
%                     BANDWIDTH             - Number of lags used in the RMK except in the case of an
%                                             infinite lag kernel.
%                     ADJUSTMENT            - De-biasing adjustment for the RK to recognize that less than
%                                             the entire sample was used to compute the 0-lag term
%                     RTPRICE               - Prices used in the RK
%                     RTTIME                - Time stamps of filtered prices used in the RK
%                     NOISEVARIANCE         - Bandi-Russell noise variance estimate. Empty if not needed to
%                                             compute the optimal bandwidth or number of points to jitter
%                     DEBIASEDNOISEVARIANCE - Bias adjusted estimate of the noise variance
%                     JITTERLAGS            - Present only if the end points are jittered, contains the
%                                             number of points that used to implement the pre-averaging.
%
% COMMENTS:
%   The only samplng scheme available is refresh time.  Hence the input SAMPLINGTYPE is not defined.
%   The only end point treatment is jittering (pre-averaging).  This is required to ensure the
%   multivariate realized kernels are positive semi-definite.
%
%   Calling realized_multivariate_kernel with a single price series:
%   [RMK] = realized_multivariate_kernel(PRICE1,TIME1,TIMETYPE,SAMPLINGINTERVAL,OPTIONS)
%
%   is the same as calling realized_kernel:
%   [RK] = realized_kernel(PRICE1,TIME1,TIMETYPE,'BusinessTime',SAMPLINGINTERVAL,OPTIONS)
%
% EXAMPLES:
%
%  See also REALIZED_COVARIANCE, REALIZED_HAYASHI_YOSHIDA, REALIZED_KERNEL, REALIZED_VARIANCE,
%  REALIZED_QUANTILE_VARIANCE, REALIZED_RANGE




k = length(varargin);
options = varargin{k};
timeType = varargin{k-2};
samplingInterval = varargin{k-1};
numPrices = (k-3)/2; % it is just how many stocks,last 3 inputs are timetype and sampling interval and option

% 0. Quick filter using refresh time to estimate the number of data points that will be available
rtPrice = realized_refresh_time_ji(timeType,varargin{1:2*numPrices}); % sychronize time of different stocks
% if samplingInterval~=1
%     rtPrice=rtPrice(1:samplingInterval:size(rtPrice,1),:);
% end
nMax = size(rtPrice(1:samplingInterval:size(rtPrice,1),:),1);

% 1. Determine optimal jitter and jitter each series [if not provided]
% and
% 2. Determine optimal number of lags for each series [if not provided]
options.filteredN = nMax;

% if isempty(options.jitterLags) && isempty(options.bandwidth)
% Need both
bandwidth = zeros(numPrices ,1);
jitterLags = zeros(numPrices ,1);
noiseVariance = zeros(numPrices ,1);
IQEstimate = zeros(numPrices ,1);
xiSquared = zeros(numPrices ,1);
for i=1:numPrices  % iterate on each stock
    % Extract price and time
    price = varargin{2*(i-1)+1};
    time = varargin{2*i};
    
    [noiseVariance(i), ~, IQEstimate(i)] = realized_noise_estimate_ji(price, time, timeType, options);
    %  bandwidth(i) = realized_kernel_bandwidth(noiseVariance(i), IQEstimate(i), options);
    % Estimation of the optimal bandwidth for use in Realized Kernels
    xiSquared(i) = noiseVariance(i)/sqrt(IQEstimate(i));
    cStar = ((12)^2/0.269)^(1/5);
    bandwidth(i) = cStar * (xiSquared(i))^(2/5) * nMax^(3/5);
    
    
    if bandwidth>(.25*nMax)
        bandwidth = round(.25*nMax);
        warning('oxfordRealized:realizedKernelLength','The estimated bandwidth requires a lag length larger than 25%% of the available data.Band width has been truncated to 25%% of data.');
    end
    kernel = options.kernel; % It's useless, because this function only works for nonflatparzen
    jitterLags(i) = realized_kernel_jitter_lag_length_ji(noiseVariance(i), IQEstimate(i), nMax);
    jitterLags(i) = max(jitterLags(i),1);
end
options.bandwidth = round(mean(bandwidth));% add the calculated bandwidth into option for later use


for i=1:numPrices
    price = varargin{2*(i-1)+1};
    time = varargin{2*(i-1)+2};
    m = size(price,1);
    p0 = mean(price(1:jitterLags(i)));
    p1 = mean(price(m-jitterLags(i)+1:m));
    t0 = time(ceil(mean(1:jitterLags(i))));
    t1 = time(floor(mean(m-jitterLags(i)+1:m)));
    price = [p0;price(jitterLags(i)+1:m-jitterLags(i));p1];
    time = [t0;time(jitterLags(i)+1:m-jitterLags(i));t1];
    varargin{2*(i-1)+1} = price;
    varargin{2*(i-1)+2} = time;
end
% 3. Put into refresh time
[rtPrice,rtTime] = realized_refresh_time_ji(timeType,varargin{1:2*numPrices});

if samplingInterval~=1
    rtPrice=rtPrice(1:samplingInterval:size(rtPrice,1),:);
end
% 4. Compute the weights for nonflatparzen kernel
H = options.bandwidth;
x = (1:H)';
x = x./(H+1);
weights = (1-6*x.^2+6*x.^3).*(x>=0 & x<=1/2) + 2*(1-x).^3.*(x>1/2 & x<1);
% 5. Compute the kernel
rtReturns = diff(log(rtPrice));
rmk = realized_multivariate_kernel_core(rtReturns,weights);

diagnostics.kernel     = options.kernel;
diagnostics.bandwidth  = options.bandwidth;
diagnostics.jitterLags = jitterLags;
diagnostics.weights    = weights;
diagnostics.rtPrice = rtPrice;
diagnostics.rtTime = rtTime;





function rmk = realized_multivariate_kernel_core(returns, weights)
% Realized Multivariate Kernel core routine that computes the value of a multivariate realized kernel
% given a set of returns and weights.
%
% USAGE:
%   [RMK] = realized_kernel_core(RETURNS,WEIGHTS)
%
% INPUTS:
%   RETURNS   - m by 1 column vector of returns
%   WEIGHTS   - H by 1 column vector of kernel weights corresponding to
%                    lags 1, 2, ..., H.  H should be much smaller than m
%
% OUTPUTS:
%   RMK        - Realized multivariate kernel value
%
% COMMENTS:
%   This is a helper function for REALIZED_MULTIVARIATE_KERNEL
%
%  See also REALIZED_MULTIVARIATE_KERNEL, REALIZED_KERNEL, REALIZED_PRICE_FILTER,
%  REALIZED_KERNEL_WEIGHTS, REALIZED_KERNEL_SELECT_LAG_LENGTH, REALIZED_COVARIANCE

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 1    Date: 7/1/2008


% Get the number of returns
[m,k]=size(returns);
% Get the size of the kernel
H=length(weights);
% Use all returns when computing autocovariances
gammaH = zeros(k,k,H);
% Loop and compute the values for gamma
for i=1:H
    returnsMinus  = returns(1:m-i,:);
    returnsPlus   = returns(i+1:m,:);
    gammaH(:,:,i) = returnsMinus' * returnsPlus;
end
% Compute the 0-lag gamma using all returns
gamma0 = returns' * returns;
% Construct the kernel
rmk = gamma0;
for j=1:length(weights)
    rmk = rmk + weights(j)*(gammaH(:,:,j)+gammaH(:,:,j)');
end



