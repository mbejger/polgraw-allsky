function [PFc,mco,bbn,Ncan,Nkall,sigpark,pm] = CoinSignificanceO1(nod,GridDir,flsum,k,veto_fraction)
% COINSIGNIFICANCEO1 False alarm probability PFc of coincidences for band
% no. k  in file  flsum and coherent search time of nod sidereal days with 
% fraction = veto_fraction of bandwidth vetoed.
% GridDir is the directory where grid.bin file is located.
%

% Bandwidth
nar = 2^2;
B = 1/nar;
dt = 1/B/2;
df = 1-2^-5;

% Parameters
N = lvcparameters(0,B,nod);                                     
% Cell size
Lc = 4;
To = N*dt;     

fid = fopen(flsum);
R1 = textscan(fid,'%s %s %f %f %f%*[^\n]',k-1);
R2 = textscan(fid,'%s %s %f %f %f %f %f %f %f %f',1);
mco = R2{5};
sss = cell2mat(R2{1});
bbn = str2double(sss(1:4));
pm = str2double(sss(6));

fpo = 10.0 + df*B*bbn;
sigpark = [R2{6} R2{7} R2{8} R2{9} R2{10}]; 

for l = 1:mco
    Nt = textscan(fid,'%f %f %f',1);
    Ncan(l,1) = Nt{1};
    Ncan(l,2) = Nt{3};
end
fclose(fid);

Nku = Ncan(:,2);
Nkall = sum(Nku);

[M,fftpad,gamrn] = freadGrid('grid.bin',GridDir);

f_min = fpo; % + 1/4/pi*B;
f_max = fpo + B; % - 1/4/pi*B;
                                      
% No. of cells
Nc = NoCells(To,f_min,f_max,gamrn);

kb = find(veto_fraction(:,1) == bbn);
if ~isempty(kb)
    Nc = round((1-veto_fraction(kb,2))*Nc/Lc^4);
else
    Nc = Nc/Lc^4;
end

[PF,PFc] = FalseAlarmCfast(2,Nc,Nku,mco);
