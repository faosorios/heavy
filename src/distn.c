/* ID: distn.c, last updated 2019/08/02, F. Osorio */

#include "distn.h"

/* static functions.. */
static double dtgamma(double, double, double, double, int);
static double ptgamma(double, double, double, double, int);
static double qtgamma(double, double, double, double, int);
static int ncomp_optimal(double);
/* ..end declarations */

static double
dtgamma(double x, double shape, double scale, double truncation, int log_pdf)
{ /* density of the right truncated gamma distribution */
  double val;
  int give_log = 1, lower_tail = 1;

  val = dgamma(x, shape, scale, give_log) - pgamma(truncation, shape, scale, lower_tail, give_log);

  return (log_pdf ? val : exp(val));
}

void
pdf_tgamma(int *n, double *y, double *x, double *shape, int *nshape, double *scale,
  int *nscale, double *truncation, int *ntrunc, int *give_log)
{ /* pdf_tgamma to be called by R */
  int i, nobs = *n, na = *nshape, nb = *nscale, nt = *ntrunc, log_pdf = *give_log;

  for (i = 0; i < nobs; i++)
    y[i] = dtgamma(x[i], shape[i % na], scale[i % nb], truncation[i % nt], log_pdf);
}

static double
ptgamma(double x, double shape, double scale, double truncation, int lower)
{ /* distribution function of the right truncated gamma */
  double val;
  int log_p = 0, lower_tail = 1;

  if (x > truncation)
    x = truncation;

  val  = pgamma(x, shape, scale, lower, log_p);
  val /= pgamma(truncation, shape, scale, lower_tail, log_p);

  return val;
}

void
cdf_tgamma(int *n, double *y, double *x, double *shape, int *nshape, double *scale,
  int *nscale, double *truncation, int *ntrunc, int *lower_tail)
{ /* cdf_tgamma to be called by R */
  int i, nobs = *n, na = *nshape, nb = *nscale, nt = *ntrunc, lower = *lower_tail;

  for (i = 0; i < nobs; i++)
    y[i] = ptgamma(x[i], shape[i % na], scale[i % nb], truncation[i % nt], lower);
}

static double
qtgamma(double p, double shape, double scale, double truncation, int lower)
{ /* quantile function of the right truncated gamma */
  double val;
  int log_p = 0, lower_tail = 1;

  if (shape == 0.0)
    return 0.0;

  p  *= pgamma(truncation, shape, scale, lower, log_p);
  val = qgamma(p, shape, scale, lower_tail, log_p);

  return val;
}

void
quantile_tgamma(int *n, double *y, double *p, double *shape, int *nshape, double *scale,
  int *nscale, double *truncation, int *ntrunc, int *lower_tail)
{ /* quantile_tgamma to be called by R */
  int i, nobs = *n, na = *nshape, nb = *nscale, nt = *ntrunc, lower = *lower_tail;

  for (i = 0; i < nobs; i++)
    y[i] = qtgamma(p[i], shape[i % na], scale[i % nb], truncation[i % nt], lower);
}

double
rng_tgamma(double shape, double rate, double truncation)
{ /* random number generation from the gamma with right truncation point,
   * i.e. TG^-(a,b,t) */
  double x;

  x = rng_tgamma_standard(shape, rate * truncation) * truncation;
  return x;
}

void rand_tgamma(int *n, double *x, double *shape, int *nshape, double *scale, int *nscale,
  double *truncation, int *ntrunc)
{ /* rand_tgamma to be called by R */
  int i, nobs = *n, na = *nshape, nb = *nscale, nt = *ntrunc;

  GetRNGstate();
  for (i = 0; i < nobs; i++)
    x[i] = rng_tgamma(shape[i % na], scale[i % nb], truncation[i % nt]);
  PutRNGstate();
}

/* random number generation of right truncated Gamma distribution using mixtures.
 * Original C code from Anne Philippe (1997). */

static int
ncomp_optimal(double b)
{ /* optimal number of components for p = 0.95 fixed */
  double ans, q = 1.644853626951;
  ans = 0.25 * R_pow_di(q * sqrt(q * q + 4. * b), 2);
  return ((int) ftrunc(ans));
}

double
rng_tgamma_standard(double a, double b)
{ /* random number generation from the gamma with right truncation point t = 1,
   * i.e. TG^-(a,b,1) */
  int n, i, j;
  double x, u, y, z, yy, zz;
  double *wl, *wlc;

  n = ncomp_optimal(b);
  wl  = (double *) Calloc(n + 2, double);
  wlc = (double *) Calloc(n + 2, double);

  wl[0] = 1.0; wlc[0] = 1.0;
  for (i = 1; i <= n; i++) {
    wl[i]  = wl[i-1] * b / (a + i);
    wlc[i] = wlc[i-1] + wl[i];
  }
  for (i = 0; i <= n; i++)
    wlc[i] = wlc[i] / wlc[n];
  y = 1.0; yy = 1.0;
  for (i = 1; i <= n; i++) {
    yy *= b / i;
    y  += yy;
  }

  repeat {
    u = unif_rand();
    j = 0;
    while (u > wlc[j])
      j += 1;
    x = rbeta(a, (double) j + 1);
    u = unif_rand();
    z = 1.0; zz = 1.0;
    for (i = 1; i <= n; i++) {
      zz *= (1 - x) * b / i;
      z  += zz;
    }
    z = exp(-b * x) * y / z;
    if (u <= z)
      break;
  }
  Free(wl); Free(wlc);
  return x;
}
