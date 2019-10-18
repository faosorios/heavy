/* ID: grubbs_fit.h, last updated 2019/08/23, F. Osorio */

#ifndef HEAVY_GRUBBS_H
#define HEAVY_GRUBBS_H

#include "base.h"
#include "matrix.h"
#include "family.h"
#include "random.h"

/* structure to hold model results */
typedef struct GRUBBS_struct {
  DIMS dm;        /* dimension data info */
  FAMILY family;  /* family data and info */
  int
    *pdims;       /* dimensions */
  double
    *y,           /* measurements */
    *settings,    /* settings */
    *mu,          /* center parameters */
    *phi,         /* dispersion parameters */
    *scale,       /* scale estimate */
    *z,           /* latent characteristic */
    *residuals,   /* residuals */
    *distances,   /* mahalanobis distances */
    *weights,     /* weights for heavy-tailed distributions */
    *acov,        /* asymptotic covariance matrix */
    *control;     /* control settings for the estimation algorithm */
  int
    maxIter,      /* maximun number of iterations */
    fixShape,     /* must estimate shape parameters? */
    ndraws;       /* independent draws for Monte Carlo integration */
  double
    tolerance;    /* convergence tolerance */
} GRUBBS_struct, *GRUBBS;

/* estimation in Grubbs' model under heavy tailed distributions */
void grubbs_fit(double *, int *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *);

#endif /* HEAVY_GRUBBS_H */
