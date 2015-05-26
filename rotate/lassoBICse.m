function [ betastar ] = lassoBICse( X, y )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


    [beta, info] = lassolar(X, y);
    [~, bestIdx] = min(info.BIC);
%     ind=beta(:,bestIdx)~=0;
%     aa = sum(ind);
%     XA = X(:,ind);
    betastar = beta(:,bestIdx);



end

