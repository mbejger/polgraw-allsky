#define C_SI 299792458

typedef enum detectors {
  none,
  GB,
  NA,
  AO,
  HO,
  PR,
  TD,
  PK,
  JB,
  G3,
  G1RAD,
  G8,
  VL,
  BO,
  MO,
  NC,
  EF,
  FB,
  H1,
  H2,
  LHO,
  LLO,
  L1,
  GEO,
  G1,
  V1,
  VIRGO,
  TAMA,
  T1
} Detectors;

extern const char* const names[];

Detectors get_detector (char *);

int get_barycenter (double, Detectors, EphemerisData *,	\
		    double *, double *, double, int);
Detectors get_detector (char *);
int get_position (Detectors, double *);
double gps2mjd (double);
double sid (double, double);

