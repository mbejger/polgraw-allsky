function [PF0,PF04cor,NF,C,pfe] = FalseAlarmCfast(Cmax,Nc,Nk,noc)
% FALSEALARMCfast  Probability PF0 that there are Cmax or more coincidences by
% chance in one or more out of Nc cells among L = length(Nk) sets of candidates
% where each set l contains Nk(l) candidates, noc is maximum number of coicidences. 
% PF04cor is the false alarm probability
% assuming that the paramater space is 4-dimensional and that we obtain
% coincidences by 2^4=16 shifts of cells by 1/2 of the size of the cell in
% all possible ways
%
% A faster version of FalseAlarmC that avoids some repeated calculations
%

[PF0,NF0,pf0,C] = FalseAlarm(Cmax,Nc,Nk,noc);
[PF0a,NF0a,pf0a,Ca] = FalseAlarm(Cmax,2^1*Nc,Nk,noc);
[PF0b,NF0b,pf0b,Cb] = FalseAlarm(Cmax,2^2*Nc,Nk,noc);
[PF0c,NF0c,pf0c,Cc] = FalseAlarm(Cmax,2^3*Nc,Nk,noc);
[PF0d,NF0d,pf0d,Cd] = FalseAlarm(Cmax,2^4*Nc,Nk,noc);

for k = 1:length(C)
    pfe(k)  = sum(C(k:end));
    pfea(k) = sum(Ca(k:end)); 
    pfeb(k) = sum(Cb(k:end));
    pfec(k) = sum(Cc(k:end));
    pfed(k) = sum(Cd(k:end));
end

PF0 = 1 - (1 - pfe).^Nc;

% dim = 4
PF04cor = 1 - ( 1 - ( 2^4*pfe - (nchoosek(4,1)*pfea+nchoosek(4,2)*pfeb+nchoosek(4,3)*pfec+nchoosek(4,4)*pfed) ...
                              - (nchoosek(4,2)*pfeb+nchoosek(4,3)*pfec+nchoosek(4,4)*pfed ) ...
                              - (nchoosek(4,3)*pfec+nchoosek(4,4)*pfed) ...
                              -  pfed) ).^Nc;
 
%disp([PF0 PF04cor])
% Expected number of n-tuple coincidences
NF = Nc*C;   %Nc*pfe;
% Test
% Lv = 0:length(Nk);
% sum(Nc*C.*Lv) - sum(Nk) must be zero
