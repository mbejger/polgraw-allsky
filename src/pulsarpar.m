function  [alfa,delta,F0,F1,F2,ang,PEPOCH] = pulsarpar(pulsar)
% PULSARPAR  Parameters of pulsar hardware injections: Pulsar 2, Pulsar 3,
%            Pulsar 4, Pulsar 8, Crab, and Vela pulsar.
% Examples: [alfa,delta,F0,F1,F2,ang,PEPOCH] = pulsarpar('2');
%           [alfa,delta,F0,F1,F2,ang,PEPOCH] = pulsarpar('Crab');

%Rad to deg
deg = 180/pi;

% GPS 2 MJD convertion
%gpsS3 = 751680013;
%mjdS3 = gps2mjd(gpsS3);       %52944
%mjdlalS3 = gps2mjdlal(gpsS3); %52944.00074287037

switch  pulsar

    case '2'

        % Ephemeris of the pulsar 2
        %Epoch of period or frequency (MJD)
        PEPOCH =  52944.00074287037; %01-Nov-2003 00:00:00.000000
        % Right ascention
        alfa =   (14 + 21/60 + 01.4799511152/3600)*15/deg;
        % Declination
        delta =  (3 + 26/60 + 38.3625755849/3600)/deg;
        % Frequency
        F0 = 287.5817865;
        % First derivative
        F1 = -6.85e-14;
        % Second derivative
        F2 =  0;
        % ho phi0 psi cos(iota)
        ang = [4.0185185179e-24 4.03 -0.221788475 cos(2.76135386361)];

    case '3'
        
        % Ephemeris of the pulsar 3
        %Epoch of period or frequency (MJD)
        PEPOCH =  52944.00074287037; %01-Nov-2003 00:00:00.000000
        % Right ascention
        alfa =   (11 + 53/60 + 29.4177660885/3600)*15/deg;
        % Declination
        delta = -(33 + 26/60 + 11.7687307074/3600)/deg;
        % Frequency
        F0 = 54.4285797; 
        % First derivative
        F1 =  -7.3e-18;
        % Second derivative
        F2 =  0;
        % ho phio psi cos(iota)
        ang = [1.62770861414e-23 5.53 0.444280306  cos(1.6515e+000)];   %1.65154960827

    case '4'
        
        % Ephemeris of the pulsar 4
        %Epoch of period or frequency (MJD)
        PEPOCH =  52944.00074287037; %01-Nov-2003 00:00:00.000000
        % Right ascention
        alfa =   (18 + 39/60 + 57.0428284445/3600)*15/deg;
        % Declination
        delta = -(12 + 27/60 + 59.8485847268/3600)/deg;
        % Frequency
        F0 = 701.5816655;
        % First derivative
        F1 =  -1.27e-08;
        % Second derivative
        F2 =  0;
        % ho phio psi cos(iota)
        ang = [4.56204779635e-23 4.83 -0.647939117 cos(1.28979208946)];
       
    case '5'
        
        % Ephemeris of the pulsar 5
        %Epoch of period or frequency (MJD)
        PEPOCH = 52944.00074287037; %01-Nov-2003 00:01:04.184000
        %PEPOCH =  52944; %01-Nov-2003 00:00:00.000000
        % Right ascention
        alfa = (20 + 10/60 + 30.3939266193/3600)*15/deg;
        %alfa =  5.281831296;
        % Declination
        delta = -(83 + 50/60 + 20.9035791210/3600)/deg;
        %delta = -1.463269033;
        % Frequency
        F0 = 26.40416218;   
        % First derivative
        F1 = -4.03e-18/2;
        % Second derivative
        F2 =  0;
        % ho phio psi cos(iota)
        ang = [4.84996486239e-24  2.23 -0.363953188 cos(1.08945639689)];
        
    case '8'
        
        % Ephemeris of the pulsar 8
        %Epoch of period or frequency (MJD)
        PEPOCH =  52944.00074287037; %01-Nov-2003 00:00:00.000000
        % Right ascention
        alfa =   (23 + 25/60 + 33.4997197871/3600)*15/deg;
        % Declination
        delta = -(33 + 25/60 + 6.6608320859/3600)/deg;
        % Frequency
        F0 = 97.15415925;  
        % First derivative
        F1 =  -4.325e-09;
        % Second derivative
        F2 =  0;
        % ho phio psi cos(iota)
        ang = [1.58762900137e-23 5.89 0.170470927 cos(1.49682623373)];
        
    case 'Vela'

        % Ephemeris of the Vela pulsar
        %PSRJ = 'PSR J0835-4510'
        %Epoch of period or frequency (MJD)
        PEPOCH =  54620; %03-Jun-2008 00:00:00.000000
        % Right ascention
        %RAJ             08:35:20.6742931
        alfa =   (8 + 35/60 + 20.6742931/3600)*15/deg;
        % Declination
        %DECJ           -45:10:36.52117
        delta = -(45 + 10/60 + 36.52117/3600)/deg;
        % Frequency
        F0 = 11.1905728445107;
        % First derivative
        F1 =  -1.55774281157927e-11;
        % Second derivative
        F2  =  4.03370032122227e-22;
        % ho phio psi cos(iota)
        ang = [1 pi/7 130.63/deg  cos(63.6/deg)];
        
        
    case 'Crab1'       
        % Ephemeris of the Crab pulsar (VSR2)
        %PSRB = 'PSR B0531+21';
        %PSRJ = 'PSR J0534+2200';
        %Epoch of period or frequency (MJD)
        PEPOCH = 40000; %16-Jan-2000 07:39:21.600000
        %Right ascension (J2000)
        alfa = (5 + 34/60 + 31.973/3600)*15/deg;
        %Declination (J2000)
        delta = (22 + 00/60 + 52.06/3600)/deg;
        F0 =  30.2254370;                %5.000e-10
        F1 =  -3.86228E-10;              %2.000e-15
        F2 =   1.2426E-20;               %1.000e-24
        %F3 =  -0.64E-30;
        % ho phio psi cos(iota)
        ang = [1 pi/7 124/deg  cos(61.3/deg)];
        
    
    case 'Crab2'       
        % Ephemeris of the Crab pulsar (VSR4)
        %PSRB = 'PSR B0531+21';
        %PSRJ = 'PSR J0534+2200';
        %Epoch of period or frequency (MJD)
        PEPOCH = 54632.00000022558;
        %Right ascension (J2000)
        alfa = (5 + 34/60 + 31.97232/3600)*15/deg;
        %Declination (J2000)
        delta = (22 + 00/60 + 52.069/3600)/deg;
        F0 =  29.74654212201602;
        F1 =  -3.719908752949545e-10;
        F2 =  1.175816039273666e-20;
        %F3 =  -0.64E-30;
        % ho phio psi cos(iota)
        ang = [1 pi/7 124/deg  cos(61.3/deg)]; 
        
    otherwise
        disp('Unknown pulsar')
end

