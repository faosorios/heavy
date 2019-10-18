/* ID: distn.h, last updated 2019/08/02, F. Osorio */

#ifndef HEAVY_DISTN_H
#define HEAVY_DISTN_H

#include "base.h"

/* dpqr-functions for the right truncated gamma distribution (to be called by R) */
extern void pdf_tgamma(int *, double *, double *, double *, int *, double *, int *, double *, int *, int *);
extern void cdf_tgamma(int *, double *, double *, double *, int *, double *, int *, double *, int *, int *);
extern void quantile_tgamma(int *, double *, double *, double *, int *, double *, int *, double *, int *, int *);
extern void rand_tgamma(int *, double *, double *, int *, double *, int *, double *, int *);

/* right truncated Gamma distribution */
extern double rng_tgamma(double, double, double);
extern double rng_tgamma_standard(double, double);

#endif /* HEAVY_DISTN_H */
