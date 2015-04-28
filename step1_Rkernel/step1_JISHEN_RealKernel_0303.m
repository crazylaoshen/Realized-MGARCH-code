clc;
clear;
% Read in data in cell array
% take two stocks for demonstration
stocknames = {'ABT','ACN','AEP','BAX'};
K = length(stocknames); % number of stocks(variance)
year = {'2006','2007','2008','2009','2010','2011','2012'};
Y = length(year);
cellData = cell(251,K); % preallocate cell for store data
RK = cell(Y,1);
% RK = zeros(K,K,Y*251); % preallocate cell for store kernel based on # of years
closePrice = cell(Y,1);

%  Note: There are 3 loops, the out loop is based on # of years, 2nd loop based on #
%  of stocks, the inner loop based on # of days of each year
for y = 1:Y
    
    
for i = 1:K  % loop over different stocks
    filename = [stocknames{i} year{y} '.csv']; % get file names for each stock
    Data = csvread(filename,1,0);  % read data
    dateind = unique(Data(:,1));  % get unique date for index
    for j = 1:length(dateind)
       
        oneDay = Data(Data(:,1) == dateind(j),[2,3]);  % extract data for one day
        [~, order] = sort(oneDay(:,1)); % sort the data based on timestamp
        oneDay = oneDay(order,:);
        
% convert the timestamps format to the one can be used in the calculation
% of kernel. HHMMSSFFF to HHMMSS.FFF  eg. 93000123 to 93000.123
         oneDay(:,1) = oneDay(:,1)/1000; 
         
         
% Eliminate the same time record, if there are several same time record, get the median one.
% The input are [price,time], check the function realized_conpute_median for details         
[medianPrice,medianTime] = realized_compute_median(oneDay(:,2),oneDay(:,1));
cellData{j,i}=[medianTime,medianPrice]; % store data in cell, each cell represent one day,j is the day,i is the stock type
n = length(medianPrice);

closePrice{y,1}(j,i) = medianPrice(n);
% closePrice(251*(y-1)+j,i) = medianPrice(n);  % check why closePrice last obs are zero 251 is suspicious



    end
end
% Now we can calculate the realized kernel with the MFE toolbox 
 for d = 1:length(dateind) % loop over every day inside a year
    options = realized_options('Multivariate Kernel');
% Note to change numbers of stocks
rmk = realized_multivariate_kernel(cellData{d,1}(:,2),cellData{d,1}(:,1) ,cellData{d,2}(:,2),cellData{d,2}(:,1) ,cellData{d,3}(:,2),cellData{d,3}(:,1) ,...
    cellData{d,4}(:,2),cellData{d,4}(:,1) ,'wall',1,options);

RK{y,1}(:,:,d) = rmk;
% RK(:,:,(y-1)*251+d) = rmk;

end

end
closePrice = cat(1,closePrice{:});
RK = cat(3,RK{:});
savefile = 'RK.mat';
save(savefile, 'RK');
savefile = 'closePrice.mat';
save(savefile, 'closePrice');