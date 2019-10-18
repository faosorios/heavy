/* ID: lm_fit.h, last updated 2019/08/05, F. Osorio */

#ifndef HEAVY_LM_FIT_H
#define HEAVY_LM_FIT_H

#include "base.h"
#include "matrix.h"
#include "family.h"
#include "random.h"

/* structure to hold model results */
typedef struct LM_struct {
  DIMS dm;        /* dimension data info */
  FAMILY family;  /* family data and info */
  int
    *pdims;       /* dimensions */
  double
    *y,           /* responses */
    *x,           /* model matrix */
    *settings,    /* settings */
    *coef,        /* coefficient estimates */
    *scale,       /* scale estimate */
    *fitted,      /* fitted values */
    *resid,       /* residuals */
    *distances,   /* mahalanobis distances */
    *weights,     /* weights for heavy-tailed distributions */
    *acov,        /* coefficients covariance matrix */
    *control;     /* control settings for the estimation algorithm */
  int
    maxIter,      /* maximun number of iterations */
    fixShape,     /* must estimate shape parameters? */
    ndraws;       /* independent draws for Monte Carlo integration */
  double
    tolerance;    /* convergence tolerance */
} LM_struct, *LM;

/* estimation in uni/multivariate linear models under heavy tailed distributions */
extern void  lm_fit(double *, double *, int *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *);
extern void mlm_fit(double *, double *, int *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *);

#endif /* HEAVY_LM_FIT_H */
