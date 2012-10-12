#include <math.h>
#include "auxi.h"

int nod, N, Nv, nfft, s, nd, fftpad, interpftpad;
double dt, B, oms, deg, c, AU, epsma;
double a, f, b, Omegar, omr, SIDday, TAIday;
double ephi, elam, eheight, egam, epsi, r;
double yr, tau_min, alfa, Smax;
double c1, c2, c3, c4, c5, c6, c7, c8, c9;

void
lvcvirgo (double fpo) { /* Offset frequency fpo */

  /* Data sampling time */
  dt = .5;
  /* Bandwidth */
  B = .5/dt;
  /* Dimensionless angular bandwidth */
  /* Ob = 2.*M_PI*B*dt; */ /* = pi */
  /* Dimensionless angular frequency */
  oms = 2.*M_PI*fpo*dt;
  /* Rad to deg */
  deg = 180./M_PI;

  /* Universal constants */
  /* speed of light */
  c = 299792.458;		/* km s^-1 */
  /* Astronomical unit */
  AU = 1.49597870691e8;		/* km */
  /* Average obliquity of the ecliptic */
  epsma = 84381.448/3600./deg;	/* 23.439 */
  /* Earth ellipsoid */
  a = 6378.140;
  f = 298.257;
  b = a*(1.-1./f);
  /* Earth angular velocity Omega_r  */
  Omegar = 7.2921151467064e-5;	/* radians/second */
  omr = Omegar*dt;
  /* Sideral day */
  SIDday = 2.*M_PI/Omegar;
  /* TAI day */
  TAIday = 86400.;

  /* Geographical location of the Virgo detector */
  /* Geographical latitude phi in radians */
  ephi = (43.+37./60.+53.0880/3600.)/deg;
  /* Geographical longitude in radians */
  elam = (10.+30./60.+16.1885/3600.)/deg;
  /* Height h above the Earth ellipsoid in meters */
  eheight = 53.238;
  /* Orientation of the detector gamma */
  /*  egam = (90. - (19.0+26./60.0))/deg; */
  egam = (135. - (19.0+25./60.0+57.96/3600.))/deg;
  /* Equatorial component r of the radius vector */
  epsi = atan(b*tan(ephi)/a);
  r = a*cos(epsi) + eheight*cos(ephi)/1000.;

  /* Observation time: no. of days */
  nod = 2;
  /* No. of data points */
  Nv = N = round (nod*SIDday/dt);

  /* length of FFT */
  nfft = 1 << (int)ceil(log(N)/log(2.));
  /* No. of spindowns */
  s = 1;
  /* 1 year [s] */
  yr = 365.25*TAIday;
  /* Minimum spindown time [s] */
  tau_min = 1000.*yr;
  /* Maximum spindown (1000 years) [angular, dimensionless] */
  Smax = 2.*M_PI*(fpo+B)*dt*dt/(2.*tau_min);

  /* False alarm probability */
  alfa = .01;
  /* Degree of freedom */
  nd = 2;		/* 2*nd = no. of degrees of freedom for Chi^2 */
  /* Zero padding */
  fftpad = 2;
  interpftpad = 2;

} /* lvcvirgo() */

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
