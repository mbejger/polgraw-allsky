%function CoinSignificanceO1All(GridDir,OutDir)
% COINSIGNIFICANCEO1ALL False alarm probability of coincidences 
% GridDir is the directory where grid.bin file is located.
% OutDir is the directory where results are sotored.
% Example:   CoinSignificanceO1All('C:\Data\010','files')
%

% Files with coincidence results from coincidence code
flsum{1} = [OutDir '/summary_0000_0242.txt'];
flsum{2} = [OutDir '/summary_0243_0328.txt'];
flsum{3} = [OutDir '/summary_0329_0447.txt'];
flsum{4} = [OutDir '/summary_0448_0520.txt'];
flsum{5} = [OutDir '/summary_0521_0599.txt'];
flsum{6} = [OutDir '/summary_0600_0699.txt'];
flsum{7} = [OutDir '/summary_0700_0759.txt'];
flsum{8} = [OutDir '/summary_0760_0784.txt'];
flsum{9} = [OutDir '/coinc-summary-global-sort_5_10.txt'];

% Bandwidth
nar = 2^2;
B = 1/nar;
%dt = 1/B/2;
df = 1-2^-5;

% No. of sidereal days
nod = 6;

% Load file with fraction of vetoed bandwidth for each band
veto = load([OutDir '/veto_fraction_frame11.txt'],'-ascii');

fpo = []; 
PFce = [];
noc = []; bbn = []; Nkall = []; sigpar = []; pm = [];

for  l = 1:length(flsum)

fid = fopen(flsum{l});
R = textscan(fid,'%s %s %f %f %f%*[^\n]');
fclose(fid);

for k = 1:length(R{1})
   disp(k)
   [PFc{k},nocl(k),bbnl(k),Ncan{k},Nkalll(k),sigpark,pml(k)] = ...
       CoinSignificanceO1(nod,GridDir,flsum{l},k,veto);
   sigparl(k,:) = sigpark;
   if isinf(PFc{k}(end))
       PFcel(k) = 1;
   else
       PFcel(k) = PFc{k}(end);
   end
end

% Frequency of lower edge of the bandwidth 
fpol = 10.0 + df*B*bbnl;

fpo = [fpo fpol];
PFce = [PFce PFcel];
noc = [noc nocl]; bbn = [bbn bbnl]; Nkall = [Nkall Nkalll];
sigpar = [sigpar; sigparl];
pm = [pm pml];

clear  fpol PFcel nocl bbnl Nkalll sigparl pml
end

% Plot coincicdence false alarm probability against frequency
figure
plot(fpo,PFce,'or')
grid on
xlabel('Frequency [Hz]')
ylabel('False alarm probability')

print('-djpeg',[OutDir '/FA_6d_07.png'])

[sortPFce,isort] = sort(PFce);

significanceO1 = [bbn(isort)' fpo(isort)' PFce(isort)' noc(isort)'  Nkall(isort)' sigpar(isort,:) pm(isort)'];

save([OutDir '/significanceO1.txt'],'significanceO1','-ascii')

% Outliers - coincidence with FAP less than 10%
ksig = find(sortPFce < 0.1);

out1 = significanceO1(ksig,:);
save ([OutDir '/outliersO1.txt'],'out1','-ascii')
