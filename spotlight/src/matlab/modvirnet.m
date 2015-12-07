function [a,b,as,bs] = modvirnet(sinal,cosal,sindel,cosdel,phir,omrt,ephi,egam)
% MODVIR  amplitude modulation functions for ground based laser
% interferometric detector
%

c1 =   sin(2*egam)*(1 + sin(ephi)^2)/4;
c2 = - cos(2*egam)*sin(ephi)/2;
c3 =   sin(2*egam)*sin(2*ephi)/2;
c4 = - cos(2*egam)*cos(ephi);
c5 = 3*sin(2*egam)*cos(ephi)^2/4;
c6 =   cos(2*egam)*sin(ephi);
c7 =   sin(2*egam)*(1 + sin(ephi)^2)/2;
c8 =   cos(2*egam)*cos(ephi);
c9 =   sin(2*egam)*sin(2*ephi)/2;

cosmodf = cos(omrt); 
sinmodf = sin(omrt);

cosalfr = cosal*cos(phir) + sinal*sin(phir);
sinalfr = sinal*cos(phir) - cosal*sin(phir);
c2d = cosdel^2;
c2sd = sindel*cosdel;

c = cosalfr*cosmodf + sinalfr*sinmodf; 
s = sinalfr*cosmodf - cosalfr*sinmodf;

c2s = 2*c.^2; cs = c.*s;

% Modulation factors
a = c1*(2 - c2d)*c2s + c2*(2 - c2d)*2*cs + c3*c2sd*c + c4*c2sd*s - ...
    c1*(2 - c2d) + c5*c2d;
b = c6*sindel*c2s + c7*sindel*2*cs + c8*cosdel*c + c9*cosdel*s - ... 
    c6*sindel;
