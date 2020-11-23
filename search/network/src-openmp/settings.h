#ifndef __SETTINGS_H__
#define __SETTINGS_H__

#include <fftw3.h>
#include "struct.h"

#define NPAR 5 		/* no. of trigger parameters */

#define INT 1		/* simplest interpolation */
#define FFT 2		/* refined (fft) interpolation */

#define NAV 4096    /* to define nmin and nmax on the edges; multiplied by B=0.5/dt */   

#define RAD_TO_DEG (180/M_PI) // = 180/pi

//constants
#define C_SPEED_OF_LIGHT 299792.458 // in km/s
#define C_AU 1.49597870691e8	// Astronomical unit, km

#define C_EPSMA (84381.448/3600./RAD_TO_DEG)
//#define C_EPSMA 0.409092804222328965124688693322241306304931640625
// Average obliquity
// of the ecliptic: 23.439

//Earth ellipsoid
#define C_ELLIPSOID_A 6378.140
#define C_ELLIPSOID_F 298.257
#define C_ELLIPSOID_B (C_ELLIPSOID_A * ( 1. - 1. / C_ELLIPSOID_F ) )
//6356.755288157528

#define C_OMEGA_R 7.2921151467064e-5
#define C_SIDDAY (2.*M_PI/C_OMEGA_R) 
// 86164.09890369719 // 2.*M_PI/Omegar      // Sideral day
#define C_TAIDAY  86400.				            // TAI day

#define C_YEARSEC (365.25*C_TAIDAY)
//31557600.0          // year in seconds = 365.25 * 86400


void search_settings( Search_settings *sett );

void detectors_settings( Search_settings *sett, Command_line_opts *opts );

void rogcvir( Detector_settings *ifo ); 

void modvir(
	    double sinal,
	    double cosal, 
	    double sindel, 
	    double cosdel,	
	    int Np, 
	    Detector_settings *ifo, 
	    Aux_arrays *aux);

int lineph( double, double *, double *, char *, char *, int );

// Lines and excluded regions treatment
void narrow_down_band( Search_settings* sett,  Command_line_opts *opts );

int read_lines( Search_settings *sett, Command_line_opts *opts );

int line_in_band( double* fl, double* fr, Search_settings* sett );

void lines_veto_fraction(Search_settings* sett, int lf, int le, int vflag);

// Coincidences 
void read_trigger_files( Search_settings *sett,
			 Command_line_opts_coinc *opts, 
			 Candidate_triggers *trig);

#endif
