function [sgnlo,pm,sgnlol,sgnloll] = JobVirSigparLin(epsm,rrn,Smax,oms,N,nospnd)
% JOBVIRSIGPARLIN  Generate signal parameters
%
% sgnlo(1) - angular frequency
% sgnlo(2) - angular spin down
% sgnlo(3) - declination
% sgnlo(4) - right ascension
% sgnlo(5:8) - amplitudes
%

%rand1 = 0.476559;

%warning off

sinal = 1i;
cosal = 1i;
num = 0;
randn('state',sum(100*clock));
rand('twister',sum(100*clock));

om1next = 2*pi;

while  ~isreal(sinal) || ~isreal(cosal) || om1next > pi*((rrn(2) - rrn(1)) + rrn(1)) 
  num = num + 1;
  disp(num);
  % Generate parameters - linear phase parametrization
  % Angular frequency
  om0 =  pi*(rand(1)*(rrn(2) - rrn(1)) + rrn(1));
  % Interference at 130Hz in 50_030 data file
  %om0 = 130.001*2*pi/2 - oms;
  
  % Frequency derivative
  if nospnd == 1
      om1 = 0;
  else
      om1 = -Smax*rand(1);
      om1next = om0 + 2*om1*N;
  end
  
  cof = oms + om0;
  be1 = (2*rand(1) - 1); be2 = (2*rand(1) - 1);
  pm = round(rand(1)) + 1;

  % Right ascension and declination
  [sinal,cosal,sindel,cosdel,alfa,delta] = ...
        lin2ast(be1,be2,pm,epsm);
    
  % Amplitude (random phase and polarization)
  A = camrnaut;
  % Signal parameters - linear phase parametrization
  sgnlol = [om0 om1 cof*be1 cof*be2 A];
  sgnloll = [om0 om1 be1 be2 A];
  % Signal parameters - astrophysical parametrization
  sgnlo = [sgnlol(1:2) delta alfa sgnlol(5:8)];
end

% Test
%[be1t,be2t,pmt] = ast2lin(alfa,delta,epsm);
%disp([be1-be1t be2-be2t pm-pmt]); must be zero
