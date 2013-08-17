#include <math.h>
#include "auxi.h"

int nod, N, Nv, nfft, s, nd, fftpad, interpftpad;
double dt, B, oms, deg, c, AU, epsma;
double a, f, b, Omegar, omr, SIDday, TAIday;
double ephi, elam, eheight, egam, epsi, r;
double yr, tau_min, alfa, Smax;
double c1, c2, c3, c4, c5, c6, c7, c8, c9;

void
settings (double fpo, char *ifo) {
     
                                // Offset frequency fpo
                                // ifo: detector

    dt = 0.5;                   // data sampling time 
    B = 0.5/dt;                 // Bandwidth
    oms = 2.*M_PI*fpo*dt;       // Dimensionless angular frequency
    deg = 180./M_PI;            // Rad to deg

    // Universal constants
  
    c = 299792.458;                 // speed of light, km/s 
    AU = 1.49597870691e8;           // Astronomical unit, km
    epsma = 84381.448/3600./deg;	// Average obliquity 
                                    // of the ecliptic: 23.439 
    // Earth ellipsoid

    a = 6378.140;
    f = 298.257;
    b = a*(1.-1./f);

    Omegar = 7.2921151467064e-5;	// Earth angular velocity Omega_r, rad/s 
    omr = Omegar*dt;
    SIDday = 2.*M_PI/Omegar;        // Sideral day
    TAIday = 86400.;                // TAI day

    nod = 2;                        // Observation time in days
    Nv = N = round (nod*SIDday/dt); // No. of data points

    nfft = 1 << (int)ceil(log(N)/log(2.));    // length of FFT
    s = 1;                          // No. of spindowns
    yr = 365.25*TAIday;             // 1 year in sec.
    tau_min = 1000.*yr;             // Minimum spindown time in sec.
    Smax = 2.*M_PI*(fpo+B)*dt*dt/(2.*tau_min); // Maximum spindown (1000 years) 
                                               // [angular, dimensionless]
     

    alfa = .01;       // False alarm probability   
    nd = 2;		// Degree of freedom, 2*nd = no. of degrees of freedom for Chi^2
    fftpad = 2;   // Zero padding 
    interpftpad = 2;

    // The detectors' settings 
 
    if(ifo[0]=='V') { 
    
        // Virgo detector 

        ephi = (43.+37./60.+53.0880/3600.)/deg;     // Geographical latitude phi in rad
        elam = (10.+30./60.+16.1885/3600.)/deg;     // Geographical longitude in rad
        eheight = 53.238;                           // Height h 
                                                    // above the Earth ellipsoid in meters
        egam = (135. - (19.0+25./60.0+57.96/3600.))/deg; // Orientation of the detector gamma

        epsi = atan(b*tan(ephi)/a);                 
        r = a*cos(epsi) + eheight*cos(ephi)/1000.;  // Equatorial component r 
                                                    // of the radius vector

    } 

} // settings 

int
rogcvir (void) {

  /* Coefficients of the amplitude modulation functions
     for VIRGO detector
  */ 

// In the notation of Phys. Rev. D 58, 063001 (1998): 
// ephi = lambda (geographical latitude phi in radians)  
// egam = gamma (orientation of the detector)
//
// (see modvir function in JobNAllSky-common.c 
// for full Eqs. 12 and 13)

  c1 = .25*sin(2.*egam)*(1+sqr(sin(ephi)));
  c2 = -.5*cos(2.*egam)*sin(ephi);
  c3 = .5*sin(2.*egam)*sin(2.*ephi);
  c4 = -cos(2.*egam)*cos(ephi);
  c5 = .75*sin(2.*egam)*sqr(cos(ephi));
  c6 = cos(2.*egam)*sin(ephi);
  c7 = .5*sin(2.*egam)*(1.+sqr(sin(ephi)));
  c8 = cos(2.*egam)*cos(ephi);
  c9 = .5*sin(2.*egam)*sin(2.*ephi);

  return 0;
} /* rogcvir() */
