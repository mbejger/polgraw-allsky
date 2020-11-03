#ifdef USE_LAL
#include <math.h>
#include <strings.h>
#include <lal/LALBarycenter.h>
#include <lal/LALInitBarycenter.h>
#include <lal/LALBarycenter.h>
#include "EphemerisDetector.h"

int get_barycenter (double gps1, Detectors site, EphemerisData *edat,	\
		    double *DetSSB, double *rDet, double dt, int len) {
     int i;
     double gps;
     BarycenterInput baryinput;
     EarthState earth;
     EmissionTime emit;

     switch (site) {
     case H1:
     case H2:
     case LHO:
	  baryinput.site.location[0] = -2161414.92636;
	  baryinput.site.location[1] = -3834695.17889;
	  baryinput.site.location[2] = 4600350.22664;
	  break;
     case L1:
     case LLO:
	  baryinput.site.location[0] = -74276.04472380;
	  baryinput.site.location[1] = -5496283.71971000;
	  baryinput.site.location[2] = 3224257.01744000;
	  break;
     case V1:
     case VIRGO:
	  baryinput.site.location[0] = 4546374.09900000;
	  baryinput.site.location[1] = 842989.69762600;
	  baryinput.site.location[2] = 4378576.96241000;
	  break;
     case none:
     default:
	  baryinput.site.location[0] = 0.;
	  baryinput.site.location[1] = 0.;
	  baryinput.site.location[2] = 0.;
	  fprintf (stderr, "No coordinates found for detector site, ");
	  fprintf (stderr, "using defaults\n");
     }
     
     /* set positions in light seconds */
     for (i=0; i<3; i++)
	  baryinput.site.location[i] /= C_SI;

     /* main loop */
     for (i=0; i<len; i++) {
	  gps = gps1+i*dt;
	  /* split time into seconds and nanoseconds */
	  baryinput.tgps.gpsSeconds = floor(gps);
	  baryinput.tgps.gpsNanoSeconds =		\
	       (gps-baryinput.tgps.gpsSeconds)*1.e9;
	  /* perform Earth barycentring */
	  XLALBarycenterEarth (&earth, &baryinput.tgps, edat);

	  /* set source information & perform barycentring */
	  /* x component */
	  baryinput.alpha = 0.;
	  baryinput.delta = 0.;
	  baryinput.dInv = 0.;
	  XLALBarycenter (&emit, &baryinput, &earth);
	  DetSSB[3*i] = emit.deltaT/dt;
	  rDet[3*i] = emit.erot/dt;

	  /* y component */
	  baryinput.alpha = M_PI_2;
	  baryinput.delta = 0.;
	  baryinput.dInv = 0.;
	  XLALBarycenter (&emit, &baryinput, &earth);
	  DetSSB[3*i+1] = emit.deltaT/dt;
	  rDet[3*i+1] = emit.erot/dt;

	  /* z component */
	  baryinput.alpha = 0.;
	  baryinput.delta = M_PI_2;
	  baryinput.dInv = 0.;
	  XLALBarycenter (&emit, &baryinput, &earth);
	  DetSSB[3*i+2] = emit.deltaT/dt;
	  rDet[3*i+2] = emit.erot/dt;
     }
     return 0;
}

const int nnames = 29;

const char* const names[] = {
     "none",
     "GB",
     "NA",
     "AO",
     "HO",
     "PR",
     "TD",
     "PK",
     "JB",
     "G3",
     "G1RAD",
     "G8",
     "VL",
     "BO",
     "MO",
     "NC",
     "EF",
     "FB",
     "H1",
     "H2",
     "LHO",
     "LLO",
     "L1",
     "GEO",
     "G1",
     "V1",
     "VIRGO",
     "TAMA",
     "T1"
};


Detectors get_detector (char *name) {
     int i;
     for (i=0; i<nnames; i++)
	  if (strcmp (name, names[i]) == 0)
	       return i;
     return 0;
}

int get_position (Detectors detector, double *position) {
     /*    Geographical locations of gravitational wave detectors
	   position[0] - Geographical latitude in radians
	   position[1] - Geographical longitude in radians
	   position[2] - Height above the Earth ellipsoid in meters
	   position[3] - Orientation of the detector in radians 
     */

     double ephilal, elamlal, eheightlal, egamlal;
     double deg = 180./M_PI, xAzi, yAzi;

     switch (detector) {
     case L1:
     case LLO:
	  /* LIGO Livinston (Livingston, Louisiana, USA) */
	  ephilal = (30+(33+46.4196/60.)/60.)/deg;
	  elamlal = - (90+(46+27.2654/60.)/60.)/deg;  
	  eheightlal = - 6.574;
	  egamlal = 242.7165/deg;
	  break;
     case LHO:
     case H1:
     case H2:
	  /* LIGO Hanford (Hanford, Washington, USA) */
	  ephilal = (46+(27+18.528/60.)/60.)/deg;
	  elamlal = - (119+(24+27.5657/60.)/60.)/deg;  
	  eheightlal = 142.554;
	  egamlal  = 170.9994/deg;
	  break;
     case VIRGO:
     case V1:
	  /* VIRGO (Cascina/Pisa, Italy) */
	  ephilal = 0.76151183984;
	  elamlal = 0.18333805213;
	  eheightlal = 51.884;
	  xAzi = 0.33916285222;
	  yAzi = 5.05155183261;
	  egamlal = M_PI_2 - 0.5 * (xAzi + yAzi);
	  break;
     default:
	  fprintf (stderr, "Detector %s: unknown position\n", names[detector]);
	  ephilal = 0.;
	  elamlal = 0.;
	  eheightlal = 0.;
	  egamlal = 0.;
     }
     position[0] = ephilal;
     position[1] = elamlal;
     position[2] = eheightlal;
     position[3] = egamlal;
     return 0;
}

#endif
