#include <math.h>
#include "auxi.h"
#include "settings.h" 

// struct that contains search parameters 
// (declared in settings.h) 
search_parameters pars; 
  
void set_search_parameters(double fpo, char *ifo) { 

	double dt = 0.5;					// data sampling time 
    double B = 0.5/dt;                 	// Bandwidth
    pars.oms = 2.*M_PI*fpo*dt;      	// Dimensionless angular frequency

    // Universal constants
    double deg = 180./M_PI;				// Rad to deg
    double c = 299792.458;              // speed of light, km/s 
    double AU = 1.49597870691e8;        // Astronomical unit, km
    pars.epsma = 84381.448/3600./deg;   // Average obliquity 
                                        // of the ecliptic: 23.439 
    // Earth ellipsoid
    double a = 6378.140;
    double f = 298.257;
    double b = a*(1.-1./f);

    double Omegar = 7.2921151467064e-5; // Earth angular velocity Omega_r, rad/s 
    pars.omr = Omegar*dt;

    double SIDday = 2.*M_PI/Omegar;		// Sideral day
    double TAIday = 86400.;				// TAI day
    int nod = 2;                        // Observation time in days
    pars.Nv = round (nod*SIDday/dt);	// No. of data points
	pars.N = pars.Nv; 

    pars.nfft = 1 << (int)ceil(log(pars.N)/log(2.));	// length of FFT
    pars.s = 1;                          	//#mb No. of spindowns
    double yr = 365.25*TAIday;          // 1 year in sec.
    double tau_min = 1000.*yr;          // Minimum spindown time in sec.
    pars.Smax = 2.*M_PI*(fpo+B)*dt*dt/(2.*tau_min);	// Maximum spindown (1000 years) 
                                               		// [angular, dimensionless]


    //alfa = .01;       // Not used: false alarm probability   
    pars.nd = 2;						// Degree of freedom, 
										// 2*nd = no. of degrees of freedom for Chi^2
	
    pars.fftpad = 2;   					//#mb default value? Zero padding 
    pars.interpftpad = 2;				// interbinning

	// Detector-related settings
	//--------------------------

	double egam, ephi; 

	// Virgo detector 
	if(ifo == "V") {  

        ephi = (43.+37./60.+53.0880/3600.)/deg;     // Geographical latitude phi in rad
        double elam = (10.+30./60.+16.1885/3600.)/deg;     // Geographical longitude in rad
        double eheight = 53.238;                           // Height h 
                		                                   // above the Earth ellipsoid in meters

        egam = (135. - (19.0+25./60.0+57.96/3600.))/deg; // Orientation of the detector gamma

        double epsi = atan(b*tan(ephi)/a);
        double r = a*cos(epsi) + eheight*cos(ephi)/1000.;  //#mb Equatorial component r 
                                          		           // of the radius vector

	} 

	// In the notation of Phys. Rev. D 58, 063001 (1998): 
	// ephi = lambda (geographical latitude phi in radians)  
	// egam = gamma (orientation of the detector)

	// (see modvir function in JobNAllSky-common.c 
	// for full Eqs. 12 and 13)

	pars.c1 = .25*sin(2.*egam)*(1+sqr(sin(ephi)));
 	pars.c2 = -.5*cos(2.*egam)*sin(ephi);
	pars.c3 = .5*sin(2.*egam)*sin(2.*ephi);
	pars.c4 = -cos(2.*egam)*cos(ephi);
	pars.c5 = .75*sin(2.*egam)*sqr(cos(ephi));
	pars.c6 = cos(2.*egam)*sin(ephi);
	pars.c7 = .5*sin(2.*egam)*(1.+sqr(sin(ephi)));
	pars.c8 = cos(2.*egam)*cos(ephi);
	pars.c9 = .5*sin(2.*egam)*sin(2.*ephi);

}

