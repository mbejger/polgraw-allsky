function [TimePerTemplate,Nf] = BenchCpuTime(dd,bbs,DataDir,CpuTime)
% BenchCpuTime   
% Example: [TPT,Nf] = BenchCpuTime('01','001','/scratch2/BenchMark_fftpad1/',2721.00));
%

% Bandwidth [Hz]
B = 1;

%Earth angular velocity Omega_r [radians/second]
Omegar = 7.2921151467064e-5; 
%Sidereal day [s]
SIDday = 2*pi/Omegar;  
% Observation time [days]
nod = 2;
% Observation time [s]
To = nod*SIDday;

% No. of spin downs
s = 1;

%TAI day
TAIday = 86400; 
%1 year [s]
yr = 365.25*TAIday;
% Minimum spindown time [s]
tau_min = 1000*yr;

% Data directory
%DataDir = 'D:\AllSkyGaussTest\';
% Grid matrix
% Mn - normalized grid matrix
[M,fftpad,gamrn,Mn] = freadGrid('grid.bin',[DataDir dd]);

% Offset frequency
df = 1-2^-5;
bbb = str2double(bbs);
fpo = 100.0 + df*bbb;

% Volume of the parameter space
[vol,volf] = volg(To,fpo,fpo+B,s,tau_min);

% Reduced grid matrix (frequency reduced)
Mnr = Mn(2:end,2:end); 
vf = abs(det(Mnr));
% No. of filters
Nf = round(volf/vf);

TimePerTemplate = CpuTime/Nf;
