function [vol,volf,vols] = volg(To,f_min,f_max,s,tau_min)
% VOLG Volume of the intrinsics parameter space.
%   Both postive and negative spindowns
%   For s = 1, negative fdot, 2 hemispheres

% Paper III Eq.(71) x (2*pi*To)^2
cof = 4*pi^2*2^(2*s + 1)*pi^(s+2)/((s+3)*factorial(s+1));
tof = To^(s+3)*(To/tau_min)^(s*(s+1)/2);

% Volume of the intrinsic parameter space
vol = cof*tof*(f_max^(s+3) - f_min^(s+3));

% Paper III Eq.(92) x (2*pi*To)^2
coff = 4*pi^2*2^(2*s)*pi^(s+1)/(factorial(s+1));
toff = To^(s+2)*(To/tau_min)^(s*(s+1)/2);

% Volume of the intrinsic parameter space
% with frequency parameter taken out
volf = coff*toff*f_max^(s+2);

% Area of the sky
vols = 4*pi^2*pi*To^2*f_max^2;
