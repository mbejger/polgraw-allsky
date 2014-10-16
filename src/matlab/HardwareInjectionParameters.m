function   [sgnlo,f] = HardwareInjectionParameters(pulsar,oms,gps,dt)
% 

Search_Freq = 2;
% Harware injections
[alfa,delta,F0,F1,F2,ang,PEPOCH] = pulsarpar(pulsar);

% GPS times
% GPS time of the first sample of the injection data
t = gps(1);
% GPS time of the epoch PEPOCH
to = mjd2gpslal(PEPOCH);

% Interpolation of ephemeris parameters to the starting time
f = F0 + F1*(t-to) + F2*(t-to)^2/2;
fdot = F1 + F2*(t-to);
fddot = F2;

% Normalized offset frequency
%oms = 2*pi*Search_Freq*f;

% Filter parameters
sgnlo(1) = 2*pi*Search_Freq*f*dt - oms;
sgnlo(2) = 1/2*2*pi*Search_Freq*fdot*dt^2;
%sgnlo(5) = 1/6*2*pi*Search_Freq*fddot*dt^3;
sgnlo(3) = delta;
sgnlo(4) = alfa;