/* ID: mv_fit.h, last updated 2019/08/05, F. Osorio */

#ifndef HEAVY_MV_FIT_H
#define HEAVY_MV_FIT_H

#include "base.h"
#include "matrix.h"
#include "family.h"
#include "random.h"

/* structure to hold model results */
typedef struct MV_struct {
  DIMS dm;        /* dimension data info */
  FAMILY family;  /* family data and info */
  int
    *pdims;       /* dimensions */
  double
    *y,           /* data matrix */
    *settings,    /* settings */
    *center,      /* position parameter estimates */
    *Scatter,     /* scatter matrix estimate */
    *distances,   /* Mahalanobis distances */
    *weights,     /* weights for heavy tailed distributions */
    *aCov,        /* coefficients covariance matrix */
    *control;     /* control settings for estimation algorithm */
  int
    maxIter,      /* maximun number of iterations */
    fixShape,     /* must estimate shape parameters? */
    ndraws;       /* independent draws for Monte Carlo integration */
  double
    tolerance;    /* convergence tolerance */
} MV_struct, *MV;

/* estimation in linear models under heavy tailed distributions */
extern void mv_fit(double *, int *, double *, double *, double *, double *, double *, double *, double *, double *);

#endif /* HEAVY_MV_FIT_H */
