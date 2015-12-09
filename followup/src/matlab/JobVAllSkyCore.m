function  [F,sgnl,nof] = JobVAllSkyCore(xDat,DetSSB,phir,epsm,oms,omr,ephi,egam,Smax,N,nfft,rn,trl,...
    interpftpad,fftinterp,fftpad,nav,M,pm,mr,nr,spndr,or,jbf,sgf,h1)
% JobVAllSkyCore
%
nof = 0;
den = 1;
% warning off
% sgnl = []; F = [];
%N = length(xDat);
Nzeros = length(find(xDat==0));
crf0 = N/(N - Nzeros);
sig2 = var(xDat);

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
            [a,b] = modvir(sinalt,cosalt,sindelt,cosdelt,phir,omrt,ephi,egam);

            nSource  = [cosalt*cosdelt sinalt*cosdelt sindelt];
            shft = nSource*DetSSB;

            % Phase modulation
            het0 =  mod(nn*M(3,1) + mm*M(4,1),M(1,1));

            phase = het0*t + oms*shft;

            cosPH = cos(phase);  sinPH = sin(phase);
            % Matched filter
            xDatma = xDat.*a.*(cosPH - i*sinPH);
            xDatmb = xDat.*b.*(cosPH - i*sinPH);

            % Resampling
            shftf = shft - shft(1);
            if  interpftpad == 0
                % Simplest
                shft_tmp = 1 + round(t - shftf);
                i0 = find(shft_tmp <= 0) ; ie = find(shft_tmp > N-1);
                shft_tmp(i0) = []; shft_tmp(ie) = [];
                xDatma = [zeros(size(i0)) xDatma(shft_tmp) zeros(size(ie))];
                xDatmb = [zeros(size(i0)) xDatmb(shft_tmp) zeros(size(ie))];
            else
                % Refined
                xDatma = interpftmy(xDatma,interpftpad*N);
                xDatma = interp1((0:interpftpad*N-1), xDatma, interpftpad*(t-shftf),'*spline', 0);
                xDatmb = interpftmy(xDatmb,interpftpad*N);
                xDatmb = interp1((0:interpftpad*N-1), xDatmb, interpftpad*(t-shftf),'*spline', 0);
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
                    phase2 = het1*t + sgnlt(2)*(t2 + 2*t.*shft);
                    cosPH = cos(phase2);  sinPH = sin(phase2);
                    xa = xDatma.*(cosPH - i*sinPH);
                    xb = xDatmb.*(cosPH - i*sinPH);
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
                        F = real(xa).^2 + imag(xa).^2 + real(xb).^2 + imag(xb).^2 ;
                        
                        % Normalize F
                        F = FStat(F,fftpad*nfft/2,nav);
                    else
                        % For white noise use
                        F = ( real(xa).^2 + imag(xa).^2 + real(xb).^2 + imag(xb).^2 )/N/sig2;
                    end
                    % Correct for zeros
                    F = F/crf0;
                             
                    F(1:rn(1)) = 0;  F(rn(2):end) = 0;
                    
%                     figure(2)
%                     plot(F)
%                     grid on
%                     pause(1)
                    
                    save Ffftpad  F
                    
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
