function [filteredPrice,filteredTime,actualTime] = realized_price_filter_ji(price,time,timeType,samplingType,samplingInterval)
% Price filtering for computing realized variances and other
% high-frequency price based returns
%
% Revised for fixed(samplingtype)or BusinessUniform(samplingtype) and wall(timetype) time of JiShen
m=size(price,1);
timeType= lower(timeType);
samplingType=lower(samplingType);

if strcmpi(samplingType,'fixed')
    time0 = min(samplingInterval);
    time1 = max(samplingInterval);
    samplingInterval = wall2unit(samplingInterval,time0,time1);
elseif  strcmpi(samplingType,'businessuniform')
    time0 = min(time);
    time1 = max(time);
end
% Since we only work with timetype'wall'
if strcmpi(timeType,'wall')
time = wall2unit(time,time0,time1);
else
error('timeType not wall, need to redefine function.')
end

switch samplingType
    case 'businessuniform'
        % Sampling interval contains the number of samples
        indices=floor(linspace(1,m,samplingInterval));
        filteredPrice=price(indices);
        filteredTime=time(indices);
        actualTime = filteredTime;
    case 'fixed'
        filteredTime = samplingInterval;
        % Compute the filtered prices
        [filteredPrice,actualTime]=fasttimefilter(price,time,filteredTime);
    otherwise
        error('Unrecogniced SAMPLINGTYPE.')
end

% Finally convert filteredTime back to the original time format
filteredTime = unit2wall(filteredTime,time0,time1);
actualTime = unit2wall(actualTime,time0,time1);




function [filteredPrice,actualTime]=fasttimefilter(price,time,filteredTime)

% Get the size of the inputs
m=size(price,1);  % m is much bigger than n
n=size(filteredTime,1);
% This uses a 2-index algorithm so that it is fast for any realistic size.
timeIndex=1;
filteredTimeIndex=1;
% Pl hold the index values.  Initializing to 1 makes the back fill easy
pl=ones(n,1);
while timeIndex<=m && filteredTimeIndex<=n
    if time(timeIndex)<=filteredTime(filteredTimeIndex)
        pl(filteredTimeIndex)=timeIndex;
        timeIndex=timeIndex+1;
    elseif filteredTimeIndex<n
        % Increment filteredTimeIndex until filteredTime>=time(timeIndex)
        while filteredTimeIndex<n && filteredTime(filteredTimeIndex)<time(timeIndex)
            filteredTimeIndex = filteredTimeIndex + 1;
            % Last price interploation for these since there is no new
            % information
            pl(filteredTimeIndex) = pl(filteredTimeIndex-1);
        end
        % Only assign if you aren't on the last one
        if filteredTimeIndex<n
            pl(filteredTimeIndex)=timeIndex;
        end
    elseif time(timeIndex)>filteredTime(n)
        % No more to assign, so break
        break
    end
end

% Clean up any trailing ones
starter = find(diff(pl)<0) + 1;
if ~isempty(starter)
    pl(starter:length(pl)) = pl(starter-1);
end
% Use pl the index the filteredPrice and filterTimeActual
filteredPrice = price(pl);
actualTime = time(pl);

