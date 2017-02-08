function [PF,NF,pf,C] = FalseAlarm(Cmax,Nc,Nk,noc)
% FALSEALARM  Probability PF that there are Cmax or more coincidences by
% chance in one or more out of Nc cells among L = length(Nk) sets of candidates
% where each set l contains Nk(l) candidates, noc is maximum number of coicidences. 
% NF is the expected number of coincidences, pf is the probability that any given cell
% cointains Cmax or more coincidences, C(k) are probablities that any given
% cell contains  k = 0,1,...,L coincidences. 
%

% Number of time slots
L = length(Nk);
% Probability that a cell is occupied by a candidate
ee = Nk/Nc;

Nmax = min([noc+1 L]);

% C(n-1) is the probability that a cell contains 
% exactly n coincidences
for n = Cmax:Nmax %L
    %disp(n)
    P = nchoosek(ee,n);
    P = prod(P,2);
    Q = flipud(nchoosek(1-ee,L-n));
    Q = prod(Q,2);
    C(n-1) = sum(P.*Q);
end

% Probability that a cell cointains Cmax or more coincidences
pf = sum(C);

% False alarm probability
% Probability that there is Cmax or more coincidences in one or more cells
PF = 1 - (1 - pf)^Nc;

% Expected number of false alarms
NF = Nc*pf;


%cof = factorial(L)/(factorial(n)*factorial(L-n));
%cof = nchoosek(L,n);
%P = perms(ee);

% NCHOOSEK Binomial coefficient or all combinations.
%     NCHOOSEK(N,K) where N and K are non-negative integers returns N!/K!(N-K)!.
%     This is the number of combinations of N things taken K at a time.
%     When a coefficient is large, a warning will be produced indicating 
%     possible inexact results. In such cases, the result is only accurate 
%     to 15 digits for double-precision inputs, or 8 digits for single-precision
%     inputs.
%  
%     NCHOOSEK(V,K) where V is a vector of length N, produces a matrix 
%     with N!/K!(N-K)! rows and K columns. Each row of the result has K of 
%     the elements in the vector V. This syntax is only practical for 
%     situations where N is less than about 15.
%

%Description

%C = nchoosek(n,k) where n and k are nonnegative integers, returns
%n!/((n–k)! k!). This is the number of combinations of n things taken k at a time.
%C = nchoosek(v,k), where v is a row vector of length n, 
%creates a matrix whose rows consist of all possible combinations 
%of the n elements of v taken k at a time. Matrix C contains n!/((n–k)! k!) rows and k columns.
% Example:
% 
% The command nchoosek(2:2:10,4) returns the even numbers from two to ten, taken four at a time:
% 
%      2     4     6     8
%      2     4     6    10
%      2     4     8    10
%      2     6     8    10
%      4     6     8    10
