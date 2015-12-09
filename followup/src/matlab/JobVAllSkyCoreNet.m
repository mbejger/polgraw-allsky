function  [F,sgnl,nof] = JobVAllSkyCoreNet(xDat,DetSSB,phir,epsm,oms,omr,ephi,egam,Smax,N,nfft,rn,trl,...
    interpftpad,fftinterp,fftpad,nav,M,pm,mr,nr,spndr,or,jbf,sgf,h1)
% JobVAllSkyCore
%
nof = 0;
den = 1;
% warning off
% sgnl = []; F = [];
%N = length(xDat);
Nzeros = length(find(xDat{1}==0));
crf01 = N/(N - Nzeros);
Nzeros = length(find(xDat{2}==0));
crf02 = N/(N - Nzeros);
Nzeros = length(find(xDat{3}==0));
crf03 = N/(N - Nzeros);

sig21 = var(xDat{1})*crf01;
sig22 = var(xDat{2})*crf02;
sig23 = var(xDat{3})*crf03;

%disp(['Data variance = ' num2str(sig2)])
if exist([jbf '.mat'],'file') == 0
    mm = mr(1); nn = nr(1); ss = spndr(1);
    save([jbf '.mat'],'mm','nn','ss');
else
    load([jbf '.mat']);

    m = mr(1):den:mm-1;
    n = nr(1):den:nr(2);
    
%     [x,y] = meshgrid(m,n);
%     figure(h1)
%     if pm == 1
%         plot(x,y,'ob','MarkerSize',5)
%     else
%         plot(x,y,'.c')
%     end
%     hold on

    mr(1) = mm; %nr(1) = nn; spndr(1) = ss - 1;
end

coft = oms; % + pi;

t = 0:N-1; t2 = t.^2; omrt = omr*t;

for mm = mr(1):den:mr(2)
    %disp(mm);  
    save([jbf '.mat'],'mm','-append')
    for nn = nr(1):den:nr(2)
        %disp(nn); 
        save([jbf '.mat'],'nn','-append')

        % Obtain linear parameters
        al1 = nn*M(3,3) + mm*M(4,3);
        al2 = nn*M(3,4) + mm*M(4,4);
        cryt = (al1^2 + al2^2)/coft^2;
        if cryt > 1
            % disp(cryt)
            % Do nothing
        else
            % Plot search progress
            %             figure(h1)
            %             subplot(2,1,1)
            %             if pm == 1
            %                 plot(mm,nn,'or')
            %             else
            %                 plot(mm,nn,'.c')
            %             end
            %             hold on
             %disp(nn);
             [sinalt,cosalt,sindelt,cosdelt] = ...
                 lin2ast(al1/coft,al2/coft,pm,epsm);
             % Calculate declination and right ascention
             sgnlt(3) = asin(sindelt);
             sgnlt(4) = mod(atan2(sinalt,cosalt),2*pi);
            
            figure(h1)
            subplot(2,1,2)
            plot(sgnlt(4)*180/pi,sgnlt(3)*180/pi,'or')
            
            % Amplitude modulation functions
            [a1,b1] = modvirnet(sinalt,cosalt,sindelt,cosdelt,phir(1),omrt,ephi(1),egam(1));
            [a2,b2] = modvirnet(sinalt,cosalt,sindelt,cosdelt,phir(2),omrt,ephi(2),egam(2));
            [a3,b3] = modvirnet(sinalt,cosalt,sindelt,cosdelt,phir(3),omrt,ephi(3),egam(3));
            
            aa1 = sum(a1.^2); bb1 = sum(b1.^2);
            aa2 = sum(a2.^2); bb2 = sum(b2.^2);
            aa3 = sum(a3.^2); bb3 = sum(b3.^2);
            
            aa = aa1/sig21 + aa2/sig22 + aa3/sig23;
            bb = bb1/sig21 + bb2/sig22 + bb3/sig23;
              
            nSource  = [cosalt*cosdelt sinalt*cosdelt sindelt];
            
            shft1 = nSource*DetSSB{1};
            shft2 = nSource*DetSSB{2};
            shft3 = nSource*DetSSB{3};
            
            % Phase modulation
            het0 =  mod(nn*M(3,1) + mm*M(4,1),M(1,1));

            phase1 = het0*t + oms*shft1;
            phase2 = het0*t + oms*shft2;
            phase3 = het0*t + oms*shft3;
            
            cosPH1 = cos(phase1);  sinPH1 = sin(phase1);
            cosPH2 = cos(phase2);  sinPH2 = sin(phase2);
            cosPH3 = cos(phase3);  sinPH3 = sin(phase3);
            
            % Matched filter
            xDatma1 = xDat{1}.*a1.*(cosPH1 - i*sinPH1);
            xDatmb1 = xDat{1}.*b1.*(cosPH1 - i*sinPH1);
            xDatma2 = xDat{2}.*a2.*(cosPH2 - i*sinPH2);
            xDatmb2 = xDat{2}.*b2.*(cosPH2 - i*sinPH2);
            xDatma3 = xDat{3}.*a3.*(cosPH3 - i*sinPH3);
            xDatmb3 = xDat{3}.*b3.*(cosPH3 - i*sinPH3);
            
            
            % Resampling
            shftf1 = shft1 - shft1(1);
            shftf2 = shft2 - shft2(1);
            shftf3 = shft3 - shft3(1);
            
            if  interpftpad == 0
                % Simplest
                shft_tmp = 1 + round(t - shftf1);
                i0 = find(shft_tmp <= 0) ; ie = find(shft_tmp > N-1);
                shft_tmp(i0) = []; shft_tmp(ie) = [];
                xDatma1 = [zeros(size(i0)) xDatma1(shft_tmp) zeros(size(ie))];
                xDatmb1 = [zeros(size(i0)) xDatmb1(shft_tmp) zeros(size(ie))];
            else
                % Refined
                xDatma1 = interpftmy(xDatma1,interpftpad*N);
                xDatma1 = interp1((0:interpftpad*N-1), xDatma1, interpftpad*(t-shftf1),'*spline', 0);
                xDatmb1 = interpftmy(xDatmb1,interpftpad*N);
                xDatmb1 = interp1((0:interpftpad*N-1), xDatmb1, interpftpad*(t-shftf1),'*spline', 0);
                xDatma2 = interpftmy(xDatma2,interpftpad*N);
                xDatma2 = interp1((0:interpftpad*N-1), xDatma2, interpftpad*(t-shftf1),'*spline', 0);
                xDatmb2 = interpftmy(xDatmb2,interpftpad*N);
                xDatmb2 = interp1((0:interpftpad*N-1), xDatmb2, interpftpad*(t-shftf1),'*spline', 0);
                xDatma3 = interpftmy(xDatma3,interpftpad*N);
                xDatma3 = interp1((0:interpftpad*N-1), xDatma3, interpftpad*(t-shftf1),'*spline', 0);
                xDatmb3 = interpftmy(xDatmb3,interpftpad*N);
                xDatmb3 = interp1((0:interpftpad*N-1), xDatmb3, interpftpad*(t-shftf1),'*spline', 0);    
            end
            
            
            for ss = spndr(1):den:spndr(end)

                het1 = mod(ss*M(2,1),M(1,1));
                sgnlt(1) = het0 + het1;
                sgnlt(2) = ss*M(2,2) + nn*M(3,2) + mm*M(4,2);

                if  sgnlt(2)>= -Smax && sgnlt(2) <= 0
                    figure(h1)
                    subplot(2,1,1)
                    %plot(ss,sgnlt(2),'or')
                    plot(2*pi*or/(nfft*fftpad),sgnlt(2)*ones(size(or)),'or')
                    %pause
                    %tic;
                    %disp(ss);
                    nof = nof+1;
                   
                    save([jbf '.mat'],'ss','-append')
                    %phase2 = het1*t + sgnlt(2)*t2;
                    phase21 = het1*t + sgnlt(2)*(t2 + 2*t.*shft1);
                    cosPH1 = cos(phase21);  sinPH1 = sin(phase21);
                    phase22 = het1*t + sgnlt(2)*(t2 + 2*t.*shft2);
                    cosPH2 = cos(phase22);  sinPH2 = sin(phase22);
                    phase23 = het1*t + sgnlt(2)*(t2 + 2*t.*shft3);
                    cosPH3 = cos(phase23);  sinPH3 = sin(phase23);
                    
                    xa1 = xDatma1.*(cosPH1 - i*sinPH1);
                    xb1 = xDatmb1.*(cosPH1 - i*sinPH1);
                    xa2 = xDatma2.*(cosPH2 - i*sinPH2);
                    xb2 = xDatmb2.*(cosPH2 - i*sinPH2);
                    xa3 = xDatma3.*(cosPH3 - i*sinPH3);
                    xb3 = xDatmb3.*(cosPH3 - i*sinPH3);
                    
                    xa = xa1/sig21 + xa2/sig22 + xa3/sig23;
                    xb = xb1/sig21 + xb2/sig22 + xb3/sig23;
                    
                    % FFT
                    if strcmp(fftinterp,'int')
                        xa = fft(xa,nfft);
                        xa(1:2:nfft-1) = xa(1:nfft/2);
                        xa(2:2:nfft-2) = (xa(1:2:nfft-3) - xa(3:2:nfft-1))/sqrt(2);
                        xb = fft(xb,nfft);
                        xb(1:2:nfft-1) = xb(1:nfft/2);
                        xb(2:2:nfft-2) = (xb(1:2:nfft-3) - xb(3:2:nfft-1))/sqrt(2);
                    else
                        xa = fft(xa,fftpad*nfft); xa(fftpad*nfft/2+1:end)=[];
                        xb = fft(xb,fftpad*nfft); xb(fftpad*nfft/2+1:end)=[];
                    end

                    % F-statistic   
                    if ~isempty(nav)               
                        F = ( real(xa).^2 + imag(xa).^2 )/aa + ( real(xb).^2 + imag(xb).^2 )/bb;
                        
                        % Normalize F
                        F = FStat(F,fftpad*nfft/2,nav);
                    else
                        % For white noise use
                        F = ( real(xa).^2 + imag(xa).^2 )/aa + ( real(xb).^2 + imag(xb).^2 )/bb;
                    end
                    % Correct for zeros
                    %F = F/crf0;
                             
                    F(1:rn(1)) = 0;  F(rn(2):end) = 0;
                    
%                   figure(2)
%                   plot(F)
%                   grid on
%                   pause(1)
                    
                    save Fnet  F
                    
                    sig = find(F > trl);
                    %toc

                    if ~isempty(sig);
                        %disp([pm nn mm ss])
                        sgnl = VirSig(sig,F,sgnlt,sgf,nfft,fftpad,pm);
                        %disp([nof length(sgnl(:,1))])
                    else
                        sgnl = [];
                    end
                    
                    
                end
            end
        end
    end
end
