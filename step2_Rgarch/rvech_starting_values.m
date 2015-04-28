function [startingvals,lls,output_parameters]=rvech_starting_values(Y,X)
% Perform a grid search to find decent starting values for SCALAR_VT_VECH(P,Q) esimtation.  If
% starting values is user supplied (and thus nonempty), does nothing.
%
% USAGE:
%   [STARTINGVALS,LLS,OUTPUT_PARAMETERS] = ...
%        scalar_vt_vech_starting_values(STARTINGVALS,DATAAUG,P,Q,T,C);
%
% INPUTS:
%   STARTINGVALS - A vector of starting values or empty to perform a grid search
%   PARAMETERS   - A vector of vech GARCH process parameters: [alpha beta]'
%   DATA         - Augmented (by m back cast values) matrix of mean zero residuals
%   P            - Positive, scalar integer representing the number of lags of the innovation process
%   Q            - Non-negative scalar integer representing the number of lags of conditional covariance
%   T            - Length of the original data
%   C            - The unconditional covariance of the data (cov(data)
%
% OUTPUTS:
%   STARTINGVALS      - A vector of starting values (p+q) by 1
%   LLS               - A vector of log likelihoods corresponding to OUTPUT_PARAMETERS
%   OUTPUT_PARAMETERS - A matrix of alternative starting values, sorted by log likelihood
%
% COMMENTS:
%   See also SCALAR_VT_VECH

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 3    Date: 9/1/2005

%This code performs a quick grid search of potential starting values
a=[.1 .15 .2 .25 .3]; % The innovation term is usually small
b=[.4 .45 .5 .6 .65]; % The smoothing parameter is usually large
% count = 1;
% params = zeros(length(a)*length(b),2);

[A,B] = meshgrid(a,b);
c=cat(2,A',B');
params1=reshape(c,[],2);
params2 = [params1(:,2),params1(:,1)];
params = [params1 ; params2];
% for i=1:length(a)
%     alpha = a(i);
%     for k = 1:length(apgpb)
%         beta = apgpb(k)-alpha;
%         if beta<=apgpb(k)/2;
%             scale = 0.5 * apgpb(k)/alpha;
%             alpha = alpha*scale;
%             beta = apgpb(k)-alpha;% Reduce alpha
%         end
%         params(count,:) = [beta alpha];
%         count = count+1;
%     end
% end
params = unique(params,'rows');
% params = params(~all(params==0,2),:);


    lls=zeros(length(params),1);
    %Evaluate the likelihood at each parameter
    for i=1:size(params,1)
        parameters=params(i,:);
        lls(i)=scalar_rvech_likelihood(parameters,Y,X);
    end
    lls(imag(lls)~=0)=nan;
    [lls,pos]=sort(lls);
    output_parameters=params(pos,:);
    startingvals=output_parameters(1,:);

