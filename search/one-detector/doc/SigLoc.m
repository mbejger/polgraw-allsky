function   [gri,nmso] = ...
    SigLoc(M,xdat,DetSSB,phir,epsm,ampmod,Fo,sgnlol,pm,oms,omr,ephi,egam)
% SigLoc Location of the signal in the parameter space
%

N = length(xdat);
N0 = length(find(xdat == 0));
sig2 = var(xdat)*(N-1)/(N-1-N0);

sgnll = sgnlol(1:4)';

Mp = M';
nmso = Mp\sgnll;

% Match factor at the corners
gri1 = floor(nmso); 
dgri = [0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1; 
        0 0 0 0 1 1 1 1 0 0 0 0 1 1 1 1;
        0 0 1 1 0 0 1 1 0 0 1 1 0 0 1 1;
        0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1];

grin = gri1*ones(1,16) + dgri;

for k = 1:16
    %disp(k);
    sgnlln = Mp*grin(:,k);
    sgnln = sgnlln;
    cof = oms;  % + sgnlol(1);  %+pi;
    [sinalt,cosalt,sindelt,cosdelt,sgnln(4),sgnln(3)] = ...
        lin2ast(sgnlln(3)/cof,sgnlln(4)/cof,pm,epsm);
    F = funvir(sgnln,xdat,DetSSB,phir,ampmod,oms,omr,ephi,egam,sig2);
    FF_match(k) = sqrt(-F/Fo); 
    %disp(['Match: ' num2str(FF_match(k))]); 
end

[FFm,im] = max(FF_match);
gri = grin(:,im);

nmso = round(nmso);