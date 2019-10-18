/* ID: specfun.c, last updated 2019/08/05, F. Osorio */

#include "specfun.h"

/* static functions.. */
static void pg_asymp(double, double, double *);
static void pg_continued_fraction(double, double, double *);
static void pg_series_expansion(double, double, double *);
/* ..end declarations */

void
cdf_gamma_derivatives(double *x, double *shape, double *scale, double *deriv)
{ /* derivatives of pgamma: to be called by R */
  double y = *x, a = *shape, b = *scale;

  pgamma_derivative(y, a, b, deriv);
}

double
pgamma_asymp(double x, double a, double scale)
{ /* Computes the regularized incomplete gamma function using a normal approximation */
  double z;

  x *= scale;
  z  =  3.0 * sqrt(a) * (R_pow(x / a, 1.0 / 3.0) + 1.0 / (9.0 * a) - 1.0);

  return pnorm(z, 0., 1., 1, 0);
}

void
pgamma_derivative(double x, double a, double scale, double *deriv)
{ /* Computes the first derivative of the incomplete gamma function.
   * Algorithm 187: Applied Statistics 31, 1982, pp. 330-335 */
  if (a <= 0.0 || scale <= 0.0)
    return;	/* should not happen */

  x *= scale;

  if (a > 600.0)
    pg_asymp(x, a, deriv);
  else if (x > 1.0 && x >= a)
    pg_continued_fraction(x, a, deriv);
  else
    pg_series_expansion(x, a, deriv);
}

double
pgamma_1st_derivative(double x, double a, double scale)
{
  double ans, *deriv;

  deriv = (double *) Calloc(3, double);
  pgamma_derivative(x, a, scale, deriv);
  ans = deriv[1];

  Free(deriv);
  return ans;
}

static void
pg_asymp(double x, double a, double *deriv)
{ /* Computes the first and second derivatives of the regularized incomplete gamma
   * function using a normal approximation */
  double u, z, z_dot;

  z      =  3.0 * sqrt(a) * (R_pow(x / a, 1.0 / 3.0) + 1.0 / (9.0 * a) - 1.0);
  z_dot  = 0.5 * (R_pow(x / a, 1.0 / 3.0) + 1.0 / (3.0 * a) - 3.0) / sqrt(a);
  u  = (5.0 * R_pow(x / a, 1.0 / 3.0) / 3.0 - 1.0 / a - 3.0) / sqrt(a);
  u += 12.0 * R_pow(a, 1.5) * z * R_pow(z_dot, 2.0);

  deriv[0] = pnorm(z, 0.0, 1.0, 1, 0);
  deriv[1] = z_dot * dnorm(z, 0.0, 1.0, 0);
  deriv[2] = -0.25 * dnorm(z, 0.0, 1.0, 0) * u;
}

static void
pg_continued_fraction(double x, double a, double *deriv)
{ /* Computes the first derivative of the regularized incomplete gamma function
   * using a continued fraction */
  double lgam, psi, psi_dot, xlog;
  double i, b, g, p, s0;
  double c1, c2, c3, c4, c5, c6, d1, d2, d3, d4, d5, d6, e1, e2, e3, e4, e5, e6;
  double f, f_dot, f_ddot, s, s_dot, s_ddot;
  const static double eps = 1.0e-06, max_it = 10000.0, scalefactor = 1.0e+30;

  /* set constants */
  xlog    = log(x);
  lgam    = lgammafn(a);
  psi     = digamma(a);
  psi_dot = trigamma(a);

  /* eval 'factor' and its derivatives */
  f      = exp(a * xlog - lgam - x);
  f_dot  = f * (xlog - psi);
  f_ddot = f_dot * f_dot / f - f * psi_dot;

  /* Use a continued fraction expansion */
  p  = a - 1.0;
  b  = x + 1.0 - p;
  c1 = 1.0;
  c2 = x;
  c3 = x + 1.0;
  c4 = x * b;
  c5 = c6 = 0.0;
  s0 = c3 / c4;
  d1 = d2 = d3 = d5 = d6 = 0.0;
  d4 = -x;
  e1 = e2 = e3 = e4 = e5 = e6 = 0.0;
  i  = 0.0;

  while (i < max_it) {
    i++;

    p--;
    b += 2;
    g = i * p;

    c5 = b * c3 + g * c1;
    c6 = b * c4 + g * c2;

    d5 = b * d3 - c3 + g * d1 + i * c1;
    d6 = b * d4 - c4 + g * d2 + i * c2;

    e5 = b * e3 + g * e1 + 2.0 * (i * d1 - d3);
    e6 = b * e4 + g * e2 + 2.0 * (i * d2 - d4);

    if (fabs(c6) > DBL_EPSILON) {
      s = c5 / c6;

      if (fabs(s - s0) <= eps * s) {
        s_dot  = (d5 - s * d6) / c6;
        s_ddot = (e5 - s * e6 - 2.0 * s_dot * d6) / c6;
        deriv[0] = 1.0 - f * s;
        deriv[1] = -f * s_dot - f_dot * s;
        deriv[2] = -f * s_ddot - 2.0 * s_dot * f_dot - f_ddot * s;
        return;
      }
      s0 = s;
    }

    c1 = c3; c2 = c4; c3 = c5; c4 = c6;
    d1 = d3; d2 = d4; d3 = d5; d4 = d6;
    e1 = e3; e2 = e4; e3 = e5; e4 = e6;

    if (fabs(c5) > scalefactor) {
      /* re-scale terms in continued fraction if terms are large */
      c1 /= scalefactor; c2 /= scalefactor;
      c3 /= scalefactor; c4 /= scalefactor;
      d1 /= scalefactor; d2 /= scalefactor;
      d3 /= scalefactor; d4 /= scalefactor;
      e1 /= scalefactor; e2 /= scalefactor;
      e3 /= scalefactor; e4 /= scalefactor;
    }
  }

  /* must not reach here */
  warning("non-convergence in pg_continued_fraction");
  s_dot  = (d5 - s * d6) / c6;
  s_ddot = (e5 - s * e6 - 2.0 * s_dot * d6) / c6;
  deriv[0] = 1.0 - f * s;
  deriv[1] = -f * s_dot - f_dot * s;
  deriv[2] = -f * s_ddot - 2.0 * s_dot * f_dot - f_ddot * s;
}

static void
pg_series_expansion(double x, double a, double *deriv)
{ /* Computes the first derivative of the regularized incomplete gamma function
   * using a series expansion */
  double lgamma_plus_1, psi, psi_plus_1, psi_dot, psi_dot_plus_1, xlog;
  double f, f_dot, f_ddot, term, term_dot, term_ddot, sum, sum_dot, sum_ddot, p, rel;
  const static double max_it = 200.0;

  /* set constants */
  xlog = log(x);
  lgamma_plus_1 = lgamma1p(a);
  psi  = digamma(a);
  psi_plus_1 = psi + 1.0 / a;
  psi_dot = trigamma(a);
  psi_dot_plus_1 = psi_dot - 1. / (a * a);

  /* eval 'factor' and its derivatives */
  f      = exp(a * xlog - lgamma_plus_1 - x);
  f_dot  = f * (xlog - psi_plus_1);
  f_ddot = f_dot * f_dot / f - f * psi_dot_plus_1;

  /* Pearson's series expansion */
  term = 1.0;
  sum  = 1.0;
  term_dot  = 0.0;
  term_ddot = 0.0;
  sum_dot   = 0.0;
  sum_ddot  = 0.0;
  p = a;

  do {
    p++;

    rel = term_dot / term;
    term_dot  = rel - 1.0 / p;
    term_ddot =  term_ddot / term - rel * rel + 1.0 / (a * a);

    term *= x / p;
    sum += term;

    term_dot *= term;
    term_ddot = term_ddot * term + term_dot * term_dot / term;
    sum_dot  += term_dot;
    sum_ddot += term_ddot;

    if (p > max_it + a) { /* convergence of the expansion is not achieved */
      warning("non-convergence in pg_series_expansion");
      break;
    }

  } while (term > sum * DBL_EPSILON);

  deriv[0] = sum * f;
  deriv[1] = f_dot * sum + f * sum_dot;
  deriv[2] = f_ddot * sum + 2.0 * f_dot * sum_dot + f * sum_ddot;
}
