function [y,rin,rout] = vol(n)
% VOL Hypervolume y of an n-dimensional sphere
%   of unit radius.
%

y = pi^(n/2)/gamma(n/2+1);
rin = 2^n/n^(n/2)/y;
rout = 2^n/y;
