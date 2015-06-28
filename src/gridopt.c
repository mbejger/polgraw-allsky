// This source file is a part of a code for narrow-banded all-sky coarse
// search for periodic GW signals. Copyright: Virgo/POLGRAW group (2010).
//  *** Unofficial version/awaiting revision ***

#define _GNU_SOURCE
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_multimin.h>
#include <getopt.h>

#include "settings.h"
#include "auxi.h"
#include "struct.h" 

#define SQ2 0.70710678118654752440  /* 1/sqrt(2) */
#define SQ3 0.57735026918962576451  /* 1/sqrt(3) */
#define SQ5 0.44721359549995793928  /* 1/sqrt(5) */
#define SQ6 0.40824829046386301637  /* 1/sqrt(6) */

#ifndef DTAPREFIX
#define DTAPREFIX ../data
#endif

void permutation_next (int *, int, int);
int invm (const double *, int N, double *);
void isoclinic (double *, double *);
void optpar (double *, double *, double, double, double, double, \
       double *, double *, double *, double *);
double optpar1 (double *, double *, double, double, double, \
    double *, double *, double *, double *, double *);


double
gridopt (double MM, double *gamrnf, double delp0, double *Mopt, \
   double *thickness) {
  double A4base[] = {1., -1., 0., 0., 0.,
         1., 0., -1., 0., 0.,
         1., 0., 0., -1., 0.,
         -.8, .2, .2, .2, .2};
  double ortho[] = {SQ2, -SQ2, 0., 0., 0.,
        SQ6, SQ6, -2.*SQ6, 0., 0.,
        .5*SQ3, .5*SQ3, .5*SQ3, -1.5*SQ3, 0.,
        .5*SQ5, .5*SQ5, .5*SQ5, .5*SQ5, -2.*SQ5};
  double v0[] = {-.4, -.2, 0., .2, .4};
  double *vert_base, *vert, *A4, R0, mu=1., nu=1., sf=0., Atr[16];
  double *x2tau, *tau2x, e0[4], e0norm=0., MMcalc=1., vce;
  int cfmax[] = {0,0,0,0};
  double vmax[] = {0.,0.,0.,0.};

  /* Calculate Voronoi cell vertices by permuting v0[] */
  {
    int i, j, k, p0[] = {0, 1, 2, 3, 4};
    vert_base = (double *) calloc (600, sizeof (double));
    for (j=0; j<5; j++)
      vert_base[j] = v0[j];
    for (i=1; i<120; i++) {
      permutation_next (p0, 5, i);
      for (j=0; j<5; j++)
  vert_base[i*5+j] = v0[p0[j]];
    }
    vert = (double *) calloc (480, sizeof (double));
    for (i=0; i<120; i++)
      for (j=0; j<4; j++)
  for (k=0; k<5; k++)
    vert[4*i+j] += ortho[5*j+k]*vert_base[5*i+k];
    free (vert_base);
  }

  /* Calculate base vectors of A_4^* in 4-dimensional space */
  {
    int i, j, k;
    A4 = (double *) calloc (16, sizeof (double));
    for (i=0; i<4; i++)
      for (j=0; j<4; j++)
  for (k=0; k<5; k++)
A4[4*i+j] += ortho[5*j+k]*A4base[5*i+k];
  }

  /* Adjust the covering radius */
  {
    int i;
    double R;
    R0 = sqrt(1-MM*MM);
    R = R0/sqrt(.4);
    for (i=0; i<16; i++)
      A4[i] *= R;
    for (i=0; i<480; i++)
      vert[i] *= R;
  }

  /* Transformations between tau and x coordinates */
  {
    gsl_vector *eval;
    gsl_matrix *evec;
    gsl_eigen_symmv_workspace *w;
    gsl_matrix_view m;
    double gaml[16];
    int i, j;

    vce = 1.;
    for (i=0; i<16; i++)
      gaml[i] = gamrnf[i];
    eval = gsl_vector_alloc (4);
    evec = gsl_matrix_alloc (4, 4);
    w = gsl_eigen_symmv_alloc (4);
    m = gsl_matrix_view_array (gaml, 4, 4);
    gsl_eigen_symmv (&m.matrix, eval, evec, w);
    gsl_eigen_symmv_free (w);
    x2tau = (double *) calloc (16, sizeof (double));
    tau2x = (double *) calloc (16, sizeof (double));
    for (i=0; i<4; i++) {
      vce *= 1./sqrt(eval->data[i]);
      for (j=0; j<4; j++)
  x2tau[4*j+i] = evec->data[4*j+i]/sqrt(eval->data[i]);
    }
    /* volume of the 4-dimensional correlation hyperellipsoid */
    vce *= .5*M_PI*M_PI*sqr(1.-MM*MM);
    gsl_matrix_free (evec);
    gsl_vector_free (eval);
    invm (x2tau, 4, tau2x);
    for (i=0; i<4; i++) {
      e0[i] = tau2x[4*i]*delp0;
      e0norm += e0[i]*e0[i];
    }
    e0norm = sqrt (e0norm);
  }

  /* Find optimal shrinking factor */
  {
    int cmax=3;
    int i, j, k, l, m, n;
    double v[4], vnorm, cphi, cphimin, mut, nut, sftmp;
    for (i=0; i<=cmax; i++)
      for (j=-cmax; j<=cmax; j++)
  for (k=-cmax; k<=cmax; k++)
    for (l=-cmax; l<=cmax; l++) {
      vnorm = 0.;
      for (m=0; m<4; m++) {
        v[m] = i*A4[m]+j*A4[4+m]+k*A4[8+m]+l*A4[12+m];
        vnorm += v[m]*v[m];
      }
      vnorm = sqrt (vnorm);
      if (vnorm > e0norm) {
        mut = e0norm/vnorm;
        cphimin = 1.;
        for (n=0; n<120; n++) {
    cphi = 0.;
    for (m=0; m<4; m++)
      cphi += v[m]*vert[4*n+m];
    cphi /= vnorm*R0;
    if (fabs (cphi) < fabs (cphimin))
      cphimin = cphi;
        }
        nut = 1+(1-mut*mut)*cphimin*cphimin/(1-cphimin*cphimin);
        sftmp = mut*nut*nut*nut;
        if (sftmp > sf) {
    sf = sftmp;
    mu = mut;
    nu = nut;
    for (m=0; m<4; m++)
      vmax[m] = v[m];
    cfmax[0] = i;
    cfmax[1] = j;
    cfmax[2] = k;
    cfmax[3] = l;
        }
      }
    }
    for (i=0; i<4; i++)
      printf ("%15.12g\t%d\n", vmax[i], cfmax[i]);
    printf ("%g\t%g\n", mu, nu);
  }

  /* Lattice satisfying constraints #1, #2 */
  {
    int i, j, is=-1;
    double iso_e[16], iso_a[16];
    double x0=.1, x1=2.5, x2, f0, f1, fc, y;
    
    for (i=3; i>=0; i--)
      if (abs(cfmax[i])==1)
  is = i;
    if (is<0)
      printf ("Key vector too long: [%3d%3d%3d%3d]\n",
        cfmax[0], cfmax[1], cfmax[2], cfmax[3]);
    if (is>0) { /* swap base vectors */
      double tmp;
      for (j=0; j<4; j++) {
  tmp = A4[4*is+j];
  A4[4*is+j] = A4[j];
  A4[j] = tmp;
      }
    }
    for (j=0; j<4; j++)
      A4[j] = vmax[j];

    isoclinic (iso_e, e0);
    isoclinic (iso_a, A4);

    f0 = optpar1 (x2tau, iso_e, mu, nu, x0, &y, iso_a, A4, Atr, Mopt);
    f1 = optpar1 (x2tau, iso_e, mu, nu, x1, &y, iso_a, A4, Atr, Mopt);
    x2 = fmod (x1-(x1-x0)/(f1-f0)*f1, 2.*M_PI);
    while (f0*f1>0. || fabs(x1-x0)>1.e-3) {
      x0 = x1;
      f0 = f1;
      x1 = x2;
      f1 = optpar1 (x2tau, iso_e, mu, nu, x1, &y, iso_a, A4, Atr, Mopt);
      x2 = fmod (x1-(x1-x0)/(f1-f0)*f1, 2.*M_PI);
    } /* while */

    /* switch to the bisection method */
    x2 = x0+.5*(x1-x0);
    fc = optpar1 (x2tau, iso_e, mu, nu, x2, &y, iso_a, A4, Atr, Mopt);
    while (x2 != x0 && x2 != x1) {
      if (fc*f0 > 0.) {
  x0 = x2;
  f0 = fc;
      } else {
  x1 = x2;
  f1 = fc;
      }
      x2 = x0+.5*(x1-x0);
      fc = optpar1 (x2tau, iso_e, mu, nu, x2, &y, iso_a, A4, Atr, Mopt);
    }

    optpar1 (x2tau, iso_e, mu, nu, x2, &y, iso_a, A4, Atr, Mopt);
    Mopt[1] = Mopt[2] = Mopt[3] = Mopt[6] = Mopt[7] = 0.;
  }

  /* Check the minimal match, calculate thickness */
  {
    int i, j, k;
    double vv[4], MMt;
    for (i=0; i<120; i++) {
      for (j=0; j<4; j++) {
  vv[j] = 0.;
  for (k=0; k<4; k++)
    vv[j] += Atr[4*j+k]*vert[4*i+k];
      }
      MMt = 0.;
      for (j=0; j<4; j++)
  for (k=0; k<4; k++)
    MMt += gamrnf[4*j+k]*vv[j]*vv[k];
      MMt = 1.-MMt;
      if (MMt < MMcalc)
  MMcalc = MMt;
    }
    *thickness = vce/fabs (det (Mopt, 4));
  }

  return MMcalc;
} /* gridopt() */

void isoclinic (double *M, double *v) {
  /* left-isoclinic rotation, transforming {1,0,0,0} onto v/|v| */
  double vnorm=0., a, b, c, d;
  int i;
  for (i=0; i<4; i++)
    vnorm += v[i]*v[i];
  vnorm = sqrt (vnorm);
  a = v[0]/vnorm;
  b = v[1]/vnorm;
  c = v[2]/vnorm;
  d = v[3]/vnorm;
  M[0] = a;
  M[1] = -b;
  M[2] = -c;
  M[3] = -d;
  M[4] = b;
  M[5] = a;
  M[6] = -d;
  M[7] = c;
  M[8] = c;
  M[9] = d;
  M[10] = a;
  M[11] = -b;
  M[12] = d;
  M[13] = -c;
  M[14] = b;
  M[15] = a;
} /* isoclinic() */

int
swap_indx (int N, int i) {
  /* Calculate permutation  index of  the i-th permutation  of 1,...N.
     Given the (i-1)th permutation, the i-th is calculated by swapping
     elements swap_indx(N,i) with swap_indx(N,i)-1.
   */
  int p, r;

  if (N == 2)
    return 1;
  p = i/N;
  r = i%N;
  if (r)
    return p&01 ? r : N-r;
  else
    return p&01 ? 1 + swap_indx (N-1, p) : swap_indx (N-1, p);
} /* swap_indx() */

void
permutation_next (int *p, int N, int i) {
  int j, t;
  j = swap_indx (N, i);
  t = p[j];
  p[j] = p[j-1];
  p[j-1] = t;
} /* permutation_next() */

void optpar (double *x2tau, double *iso_e, double mu, double nu,  \
       double x, double y, double *iso_a, double *A4,   \
       double *Atr, double *Mopt) {
  int i, j, k;
  double Atmp[16], ar[16], sx, cx, sy, cy;

  sincos (x, &sx, &cx);
  sincos (y, &sy, &cy);

  for (i=0; i<16; i++)
    ar[i] = 0.;
  ar[0] = 1.;
  ar[5] = cy*cx;
  ar[6] = cx*sy;
  ar[7] = sx;
  ar[9] = -sy;
  ar[10] = cy;
  ar[13] = -cy*sx;
  ar[14] = -sy*sx;
  ar[15] = cx;

  for (i=0; i<4; i++)
    for (j=0; j<4; j++) {
      Atr[4*i+j] = 0.;
      for (k=0; k<4; k++)
  Atr[4*i+j] += x2tau[4*i+k]*iso_e[4*k+j];
    }
  for (i=0; i<4; i++) {
    Atr[4*i] *= mu;
    for (j=1; j<4; j++)
      Atr[4*i+j] *= nu;
  }

  for (i=0; i<4; i++)
    for (j=0; j<4; j++) {
      Atmp[4*i+j] = 0.;
      for (k=0; k<4; k++)
  Atmp[4*i+j] += Atr[4*i+k]*ar[4*k+j];
    }

  for (i=0; i<4; i++)
    for (j=0; j<4; j++) {
      Atr[4*i+j] = 0.;
      for (k=0; k<4; k++)
  Atr[4*i+j] += Atmp[4*i+k]*iso_a[4*j+k];
    }
  for (i=0; i<4; i++)
    for (j=0; j<4; j++) {
      Mopt[4*j+i] = 0.;
      for (k=0; k<4; k++)
  Mopt[4*j+i] += Atr[4*i+k]*A4[4*j+k];
    }
} /* optpar() */

double optpar1 (double *x2tau, double *iso_e, double mu, double nu, \
    double x, double *y, double *iso_a, double *A4,   \
    double *Atr, double *Mopt) {
  double f0, f1, fc;
  double y0=1., y1=2., y2;

  optpar (x2tau, iso_e, mu, nu, x, y0, iso_a, A4, Atr, Mopt);
  f0 = Mopt[7];
  optpar (x2tau, iso_e, mu, nu, x, y1, iso_a, A4, Atr, Mopt);
  f1 = Mopt[7];

  y2 = fmod (y1-(y1-y0)/(f1-f0)*f1, 2.*M_PI);
  while (f0*f1>0. || fabs (y1-y0)>1.e-3) {
    y0 = y1;
    f0 = f1;
    y1 = y2;
    optpar (x2tau, iso_e, mu, nu, x, y1, iso_a, A4, Atr, Mopt);
    f1 = Mopt[7];
    y2 = fmod (y1-(y1-y0)/(f1-f0)*f1, 2.*M_PI);
  } /* while */

  /* switch to the bisection method */
  y2 = y0+.5*(y1-y0);
  optpar (x2tau, iso_e, mu, nu, x, y2, iso_a, A4, Atr, Mopt);
  fc = Mopt[7];
  while (y2 != y0 && y2 != y1) {
    if (fc*f0>0.) {
      y0 = y2;
      f0 = fc;
    } else {
      y1 = y2;
      f1 = fc;
    }
    y2 = y0+.5*(y1-y0);
    optpar (x2tau, iso_e, mu, nu, x, y2, iso_a, A4, Atr, Mopt);
    fc = Mopt[7];
  }

  *y = y2;
  return Mopt[6];
} /* optpar1() */

int
fishmg (int N, double *r1, double *r2, double *gam, double *gamr) {
  int i, j;
  double t, tt, r1n, r2n, t1m, t2m, t3m, t4m;
  double m1=0., m2=0., m1t=0., m2t=0., m1t2=0., m2t2=0.;
  double m1m1=0., m1m2=0.,m2m2=0., w;
  const double ww[] = {1./3., 31./24., 5./6., 25./24.};

  for (i=0; i<N; i++) {
    r1n = r1[i]/(double)N;
    r2n = r2[i]/(double)N;
    t = (double)i/(double)(N-1);
    tt = t*t;
    if (i<4 || i>N-5) {
      j = i<N-1-i ? i : N-1-i;
      w = ww[j];
      m1 += w*r1n;
      m1m1 += w*r1n*r1n;
      m2 += w*r2n;
      m2m2 += w*r2n*r2n;
      m1m2 += w*r1n*r2n;
      m1t += w*r1n*t;
      m2t += w*r2n*t;
      m1t2 += w*r1n*tt;
      m2t2 += w*r2n*tt;
    } else {
      m1 += r1n;
      m1m1 += r1n*r1n;
      m2 += r2n;
      m2m2 += r2n*r2n;
      m1m2 += r1n*r2n;
      m1t += r1n*t;
      m2t += r2n*t;
      m1t2 += r1n*tt;
      m2t2 += r2n*tt;
    }
  }

  t1m = .5;
  t2m = 1./3.;
  t3m = .25;
  t4m = .2;
  m1 /= (double)(N-1);
  m2 /= (double)(N-1);
  m1t /= (double)(N-1);
  m2t /= (double)(N-1);
  m1t2 /= (double)(N-1);
  m2t2 /= (double)(N-1);
  m1m1 /= (double)(N-1);
  m1m2 /= (double)(N-1);
  m2m2 /= (double)(N-1);

  gam[0] = 1.;
  gam[5] = gam[1] = t1m;
  gam[10] = gam[2] = t2m;
  gam[15] = gam[3] = m1;
  gam[20] = gam[4] = m2;
  gam[6] = t2m;
  gam[11] = gam[7] = t3m;
  gam[16] = gam[8] = m1t;
  gam[21] = gam[9] = m2t;
  gam[12] = t4m;
  gam[17] = gam[13] = m1t2;
  gam[22] = gam[14] = m2t2;
  gam[18] = m1m1;
  gam[23] = gam[19] = m1m2;
  gam[24] = m2m2;

  for (i=0; i<4; i++)
    for (j=0; j<4; j++)
      gamr[4*i+j] = gam[5*i+j+6]-gam[5*i+5]*gam[j+1];
  return 0;
} /* fishmg() */

int main (int argc, char *argv[]) {
  int i, j, l, m, c;
  double gamrnf[16], delp0, Mopt[16], M[16], Mn[16], MM, thickness, \
    epsm, Tn[16], Tin[16], gamrn[16], covr[16], covrnf[16], Tf[16], \
    T[16], gamr[16], minimal_match;
  char dtaprefix[512], ident[16];

  Search_settings sett;

  strcpy (dtaprefix, TOSTR(DTAPREFIX));
  strcpy (ident, "");

  // default value of the minimal match 
  minimal_match = sqrt(.75) ; 

  while (1) {
    static struct option long_options[] = {
      {"ident", required_argument, 0, 'i'},
      {"data", required_argument, 0, 'd'},
      // minimal match 
      {"minimal match", required_argument, 0, 'm'},
      {0, 0, 0, 0}
    };

    int option_index = 0;
    c = getopt_long (argc, argv, "i:d:m:", long_options,  \
         &option_index);
    if (c == -1)
      break;
    switch (c) {
    case 'i':
      sprintf (ident, "_%s", optarg);
      break;
    case 'd':
      strcpy (dtaprefix, optarg);
      break;
    case 'm':
      minimal_match = atof(optarg);
      break;
    case '?':
      break;
    default:
      abort ();
    } /* switch c */
  } /* while 1 */

  search_settings(&sett);
  delp0 = 2.*M_PI*(double)sett.N/(double)sett.nfft/(double)sett.fftpad;

  {
    FILE *data;
    double x;
    char filename[512];
    sprintf (filename, "%s/DetSSB%s.bin", dtaprefix, ident);
    if ((data = fopen (filename, "r")) != NULL) {
      for (i=0; i<3*(sett.N)+1; i++)
        fread ((void *)(&x), sizeof (double), 1, data);

      fread ((void *)(&epsm), sizeof (double), 1, data);    
      fclose (data);
    } else {
      perror (filename);
      return 0;
    }
  }

  {
    double *r1, *r2, gamn[25];
    r1 = (double *) calloc (sett.N, sizeof (double));
    r2 = (double *) calloc (sett.N, sizeof (double));
    lineph (epsm, r1, r2, dtaprefix, ident, sett.N);
    fishmg (sett.N, r1, r2, gamn, gamrn);
    for (i=0; i<16; i++) {
      Tf[i] = 0.;
      Tn[i] = 0.;
      T[i] = 0.;
      gamrnf[i] = 0.;
    }
    Tf[0] = Tf[5] = 1.;
    Tf[10] = Tf[15] = 1000.;
    Tn[0] = Tn[10] = Tn[15] = 1/(double)sett.N;
    Tn[5] = 1/sqr((double)sett.N);
    for (i=0; i<4; i++)
      for (j=0; j<4; j++)
  for (l=0; l<4; l++)
    T[4*i+j] += Tf[4*i+l]*Tn[4*l+j];
    for (i=0; i<4; i++)
      for (j=0; j<4; j++)
  for (l=0; l<4; l++)
    for (m=0; m<4; m++)
      gamrnf[4*i+j] += Tf[4*i+l]*gamrn[4*l+m]*Tf[4*m+j];
  }

  MM = gridopt (minimal_match, gamrnf, delp0, Mopt, &thickness);

  invm (Tn, 4, Tin);
  invm (gamrnf, 4, covrnf);
  for (i=0; i<4; i++)
    for (j=0; j<4; j++) {
      M[4*i+j] = 0.;
      gamr[4*i+j] = 0.;
      covr[4*i+j] = 0.;
      Mn[4*i+j] = 0.;
      for (l=0; l<4; l++) {
  M[4*i+j] += Mopt[4*i+l]*T[4*l+j];
  for (m=0; m<4; m++) {
    gamr[4*i+j] += Tin[4*i+l]*gamrn[4*l+m]*Tin[4*m+j];
    covr[4*i+j] += T[4*i+l]*covrnf[4*l+m]*T[4*m+j];
  }
  Mn[4*i+j] += Mopt[4*i+l]*Tf[4*l+j];
      }
    }
  
  for (i=0; i<4; i++) {
    for (j=0; j<4; j++)
      printf ("%14g\t", M[4*i+j]);
    printf ("\n");
  }

  printf ("MM verified to be = %g\n",sqrt(MM));
  printf ("thickness = %g\n", thickness);

  {
    FILE *data;
    char filename[512];
    sprintf (filename, "%s/grid%s.bin", dtaprefix, ident);
    data = fopen (filename, "w");
    fwrite ((void *)&sett.fftpad, sizeof (int), 1, data);
    fwrite ((void *)M, sizeof (double), 16, data);
    fwrite ((void *)gamrn, sizeof (double), 16, data);
    fwrite ((void *)Mn, sizeof (double), 16, data);
    fclose (data);

  }

  return 0;
}
