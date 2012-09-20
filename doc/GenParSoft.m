function    sgnlo = GenParSoft(dd,bbs)

%dd = '50';
%bbs = '148';

DataDir = ['I:\VSR1\' dd '\'];

df = 1-2^-5;
bbb = str2double(bbs);
fpo = 100.0 + df*bbb;

%smn = round(5 + 0.1*(fpo - 100));

[oms omr ephi egam N dt nd Smax nfft epsma]  = lvcvirgo(fpo);

cd(DataDir)

% Grid
M = freadGrid('grid.bin',DataDir);
[spndr,nr,mr] = GridRange(M,oms,Smax);

cryt = 100;
while cryt > 1

pm = round(rand(1)) + 1;

spnd = round((spndr(2) - spndr(1))*rand(1));
spnd = spndr(1) + smn + spnd;
nn = round((nr(2) - nr(1))*rand(1));
nn = nr(1) + smn + nn;
mm = round((mr(2) - mr(1))*rand(1));
mm = mr(1) + smn + mm;

[spnd nn mm];

al1 = mm*M(3,3) + nn*M(4,3);
al2 = mm*M(3,4) + nn*M(4,4);

be1 = al1/oms;
be2 = al2/oms;

cryt = be1^2 + be2^2

end

[sinal,cosal,sindel,cosdel,alfa,delta] = ...
         lin2ast(be1,be2,pm,epsma);
     
% Band edges
edg = 2^12;
rn = [2*edg nfft-2*edg+1];
rrn(1) = rn(1)/nfft;  rrn(2) = rn(2)/nfft;
        
sgnlo(1) =  pi*(rand(1)*(rrn(2) - rrn(1)) + rrn(1));
sgnlo(2) = -(spnd*M(2,2) + mm*M(3,2) + nn*M(4,2));
sgnlo(3) = delta;
sgnlo(4) = alfa;

% Random phase and polarization of the signal
ph_o = 2*pi*rand(1);    %ph_o = pi/7;
psik = 2*pi*rand(1);    %psik = pi/3.3;
cin =  2*rand(1) - 1;   %cin = 0.73;
hop = (1 + cin^2)/2; hoc = cin;

iota = acos(cin);

Th = [iota ph_o psik hop hoc];

sgnlo(5) =  cos(2*psik)*hop*cos(ph_o)-sin(2*psik)*hoc*sin(ph_o);
sgnlo(6) =  sin(2*psik)*hop*cos(ph_o)+cos(2*psik)*hoc*sin(ph_o);
sgnlo(7) = -cos(2*psik)*hop*sin(ph_o)-sin(2*psik)*hoc*cos(ph_o);
sgnlo(8) = -sin(2*psik)*hop*sin(ph_o)+cos(2*psik)*hoc*cos(ph_o);

disp(sgnlo)



