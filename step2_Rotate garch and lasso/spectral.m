function [ root,Av,v ] = spectral( A )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
[v,d] = eig(A);
dd = diag(d);
ddd = sqrt(dd);
Av = diag(ddd);
root = v*Av*v';

end

