function  Nc = NoCells(To,f_min,f_max,gamrn)
% NOCELLS  No. of cells for observation time To,
% frequency range form f_min to f_max, 1 spin down, 
% and mormalized Fisher matrix gamrn.
%

%vl = volg(To,f_min,f_max,s,tau_min);
% LVC search plans (s = 1)
fdotmin = 10^-8/2;
fdotmax = 10^-9/2;
vl =  16/3*pi^5*To^5*(fdotmin+fdotmax)*(f_max^3-f_min^3);

vc = vol(4)/sqrt(det(gamrn));

% Only negative spindowns
%vl = vl/2;
% No. of cells (one hemisphere)
Nc = round(vl/vc);
