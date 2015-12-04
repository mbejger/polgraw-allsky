#include <math.h>

double gps2mjd (double tgps) {
  /* Conversion from gps time to mjd */
  /*
    t = gps2mjd(tgps)

    tgps   gps time (in seconds)
    t      mjd

    GPS time is 0 at 6-Jan-1980 0:00:00 and it is linked to TAI.
    It is offset from UTC for the insertion of the leap seconds.

    *** TO BE UPGRADED EVERY SIX MONTHS

    Leap seconds

 1981 JUL  1 =JD 2444786.5  TAI-UTC=  20.0       S + (MJD - 41317.) X 0.0      S
 1982 JUL  1 =JD 2445151.5  TAI-UTC=  21.0       S + (MJD - 41317.) X 0.0      S
 1983 JUL  1 =JD 2445516.5  TAI-UTC=  22.0       S + (MJD - 41317.) X 0.0      S
 1985 JUL  1 =JD 2446247.5  TAI-UTC=  23.0       S + (MJD - 41317.) X 0.0      S
 1988 JAN  1 =JD 2447161.5  TAI-UTC=  24.0       S + (MJD - 41317.) X 0.0      S
 1990 JAN  1 =JD 2447892.5  TAI-UTC=  25.0       S + (MJD - 41317.) X 0.0      S
 1991 JAN  1 =JD 2448257.5  TAI-UTC=  26.0       S + (MJD - 41317.) X 0.0      S
 1992 JUL  1 =JD 2448804.5  TAI-UTC=  27.0       S + (MJD - 41317.) X 0.0      S
 1993 JUL  1 =JD 2449169.5  TAI-UTC=  28.0       S + (MJD - 41317.) X 0.0      S
 1994 JUL  1 =JD 2449534.5  TAI-UTC=  29.0       S + (MJD - 41317.) X 0.0      S
 1996 JAN  1 =JD 2450083.5  TAI-UTC=  30.0       S + (MJD - 41317.) X 0.0      S
 1997 JUL  1 =JD 2450630.5  TAI-UTC=  31.0       S + (MJD - 41317.) X 0.0      S
 1999 JAN  1 =JD 2451179.5  TAI-UTC=  32.0       S + (MJD - 41317.) X 0.0      S
 2006 JAN  1 =JD 2453736.5  TAI-UTC=  33.0       S + (MJD - 41317.) X 0.0      S
 2009 JAN  1 =JD 2454832.5  TAI-UTC=  34.0       S + (MJD - 41317.) X 0.0      S
 2012 JUL  1 =JD 2456109.5  TAI-UTC=  35.0       S + (MJD - 41317.) X 0.0      S

  To update, add the mjd value of the new leapseconds, computed with
        t=s2mjd('1-jul-1981')
  from the prompt of matlab (with the real date); update also nleap.
 Version 1.0 - October 1999
 Part of Snag toolbox - Signal and Noise for Gravitational Antennas
 Copyright (C) 1999  Sergio Frasca - sergio.frasca@roma1.infn.it
 Department of Physics - Universita` "La Sapienza" - Rome
  */

  double leaptimes[] = {
    44786,
    45151,
    45516,
    46247,
    47161,
    47892,
    48257,
    48804,
    49169,
    49534,
    50083,
    50630,
    51179,
    53736,
    54832,
    56109,
    99999,
    /* the above line added by K.Borkowski to make the procedure */
    /* work for dates after the last leap second */
    0,
    0,
    0
  };
  int i, nleap = 17;
  double t0 = 44244.;
  double t;

  t = tgps/86400.+t0;
  for (i=0; i<nleap; i++) {
    if (t<leaptimes[i] + (i+1)/86400.) {
      /* This condition means that a leap second itself will be
	 converted to the same numerical value of MJD (or UTC) time as
	 the second directly following it. That is, it will belong to
	 1 January; e.g., 23:59:60.5 UTC will become
	 00:00:00.5. Instead, if we used this form: if (t
	 .lt. leaptimes[i] + i/86400d0) then it would mean
	 ascribing the leap second to the previous UTC second;
	 e.g. 23:59:60.5 UTC would become 23:59:59.5
	 UTC. Interestingly, both these two conditions would pass the
	 test of conversion using mjd2gps immediately followed by
	 reverse (gps2mjd) conversion! That's because formally the GPS
	 time does not contain this leap second (the mjd2gps, using
	 decimal MJD, cannot distinguish between 23:59:60.5 and
	 00:00:00.5 on the following day). However, converting the
	 leap second with gps2mjd and back with mjd2gps would
	 necessaily fail with both conditions
      */
      t -= i/86400.;
      break;
    }
  }
  return t;
}

double sid (double mjd, double elam) {
  /*
    SID True local sidereal time (in radians!)
    mjd  - Modified Julian date (UTC based)
    elam - Geographical longitude (rad)
    jd2000 = 2451545.0; mjd = jd - 2400000.5;
  */
  double t, phir;

  t = mjd + 2400000.5 - 2451545.;

  phir = M_PI*fmod(1.55811454652 + fmod(t+t, 2.)			\
		   + t*(.547581870159e-2 + t*(1.61549e-15-t*1.473e-24))
		   + elam/M_PI, 2.);
  if (phir < 0)
    phir += 2*M_PI;
  return phir;
}
