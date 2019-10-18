/* ID: specfun.h, last updated 2019/08/05, F. Osorio */

#ifndef HEAVY_SPECFUN_H
#define HEAVY_SPECFUN_H

#include "base.h"

/* dpqr-functions for the right truncated gamma distribution (to be called by R) */
extern void cdf_gamma_derivatives(double *, double *, double *, double *);

/* normal approximation of the incomplete gamma integral */
extern double pgamma_asymp(double, double, double);

/* derivative of the incomplete gamma integral */
extern void pgamma_derivative(double, double, double, double *);
extern double pgamma_1st_derivative(double, double, double);

#endif /* HEAVY_SPECFUN_H */
