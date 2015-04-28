function rv=realized_variance_ji(price,time,timeType,samplingType,samplingInterval)
% Estimates Realized Variance and subsampled Realized Variance
%
% USAGE:
%   [RV,RVSS] = realized_variance(PRICE,TIME,TIMETYPE,SAMPLINGTYPE,SAMPLINGINTERVAL,SUBSAMPLES)
%
% INPUTS:
%   PRICE            - m by 1 vector of high frequency prices
%   TIME             - m by 1 vector of times where TIME(i) corresponds to PRICE(i)
%   TIMETYPE         - String describing the way times are measured
%                       'wall'    24-hour clock of the form HHMMSS, e.g. 101543 or 153217
%                       'seconds' Time measured in seconds past midnight on the first day.
%                       'unit'  Unit normalized date format, e.g. .1, .234, .9
%                         Unit normalized times are more general than the other types and can be
%                         applied to data from more than one calendar day
%   SAMPLINGTYPE     - String describing the type of sampling to use when
%                        filtering PRICE
%                        'CalendarTime' - Sample in calendar time using observations separated by 
%                          SAMPLINGINTERVAL seconds
%                        'CalendarUniform' - Sample in calendar time using SAMPLINGINTERVAL
%                          observations spread uniformly between TIME(1) and TIME(m)
%                        'BusinessTime' - Sample in business (tick) time using observation separated
%                          by SAMPLINGINTERVAL ticks
%                        'BusinessUniform' - Sample in business (tick) time using observations
%                          uniformly spaced in business time.
%                        'Fixed' - Sample at specific points in time. When using fixed,
%                          SAMPLINGINTERVAL must be a n by 1 vector of times with the same TIMETYPE
%                          as TIME (i.e. seconds if TIME is in seconds)
%   SAMPLINGINTERVAL  - Scalar integer or n by 1 vector whose meaning depends on SAMPLINGTYPE
%   SUBSAMPLES        - [OPTIONAL] Scalar integer indicating the number of subsample realized
%                         variance estimators to average with the original realized variance.
%                         Subsample realized variances are based on prices uniformly spaced between
%                         the times (Calendar sampling) or ticks (Business sampling).  SUBSAMPLES should 
%                         be greater than 1, since 1 uses a single subsample and so is equivalent 
%                         to the original RV will compute a subsample realized variance using the 
%                         mid-point of the price sample points, 2 will use the original sample and 
%                         the midpoint between the original sampling times, and so on.  
%                         In general this number should be small so the subsample estimators will be "sparse".
%
% OUTPUTS:
%   RV                - Realized variance estimate
%   RVSS              - Realized variance including subsample averaging
%
% COMMENTS:
%
% EXAMPLE:
%  % 5-minute realized variance using wall time returns and fixed increments
%  fixedInterval = seconds2wall(wall2seconds(93000):300:wall2seconds(160000));
%  RV = realized_variance(PRICE,TIME,'wall','Fixed',fixedInterval)
%
%  % 10-tick realized variance
%  RV = realized_variance(PRICE,TIME,'wall','BusinessTime',10)
%
%  % 5-minute realized variance with subsampling every minute
%  fixedInterval = seconds2wall(wall2seconds(93000):300:wall2seconds(160000));
%  [RV,RVSS] = realized_variance(PRICE,TIME,'wall','Fixed',fixedInterval,5)
%
%  See also REALIZED_KERNEL, REALIZED_QUANTILE_VARIANCE, REALIZED_RANGE, REALIZED_KERNEL_CORE,
%  REALIZED_PRICE_FILTER, REALIZED_KERNEL_WEIGHTS, REALIZED_KERNEL_SELECT_LAG_LENGTH

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 1    Date: 5/1/2008

% revised edition without subsample


logPrice =log(price);
% Filter prices and compute the RV
filteredLogPrice = realized_price_filter_ji(logPrice,time,timeType,samplingType,samplingInterval);
returns = diff(filteredLogPrice);
rv=returns' * returns;


% subsampledLogPrices = realized_subsample(logPrice,time,timeType,samplingType,samplingInterval,subsamples);
% rvs = zeros(subsamples,1);
% totalCount  = 0;
% for i=1:subsamples
%     returns = diff(subsampledLogPrices{i});
%     rvs(i) = (returns' * returns);
%     if i==1
%         baseCount = length(returns);
%     end
%     totalCount = totalCount + length(returns);
% end
% rvSS = sum(rvs) * (baseCount / totalCount);
