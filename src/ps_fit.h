/* ID: ps_fit.h, last updated 2019/08/05, F. Osorio */

#ifndef HEAVY_PS_FIT_H
#define HEAVY_PS_FIT_H

#include "base.h"
#include "family.h"
#include "matrix.h"
#include "optim.h"

/* structure to hold model results */
typedef struct SPLINE_struct {
  DIMS dm;        /* dimension data info */
  FAMILY family;  /* family data and info */
  int
    *pdims;       /* dimensions */
  double
    *y,           /* responses */
    *b,           /* model matrix */
    *half,        /* half of penalty matrix */
    *settings,    /* settings */
    *coef,        /* coefficient estimates */
    *scale,       /* scale estimate */
    *lambda,      /* smoothing parameter */
    *edf,         /* effective degrees of freedom */
    *GCV,         /* weighted Generalized Cross-Validation */
    *pen,         /* roughness penalty */
    *logLik,      /* penalized log-likelihood */
    *fitted,      /* fitted values */
    *resid,       /* residuals */
    *distances,   /* mahalanobis distances */
    *weights,     /* weights for heavy-tailed distributions */
    *control;     /* control settings for the estimation algorithm */
  int
    maxIter,      /* maximun number of iterations */
    fixShape,     /* must estimate shape parameters? */
    ndraws,       /* independent draws for Monte Carlo integration */
    ncycles;      /* number of cycles for the nested EM algorithm */
  double
    tolerance;    /* convergence tolerance */
} SPLINE_struct, *SPLINE;

/* GCV info required by the minimizer */
typedef struct GCV_pars {
  DIMS dm;
  double edf, WGCV, pen, RSS;
  double *u, *r, *s, *z, *rhs, *a, *fitted, *resid;
} GCV_pars, *GCVpars;

/* routines for P-spline estimation */
extern void ps_fit(double *, double *, double *, int *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *);
extern void ps_combined(double *, double *, double *, int *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *);

#endif /* HEAVY_PS_FIT_H */
