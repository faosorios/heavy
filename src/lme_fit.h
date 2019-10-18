/* ID: lme_fit.h, last updated 2019/08/05, F. Osorio */

#ifndef HEAVY_LME_FIT_H
#define HEAVY_LME_FIT_H

#include "base.h"
#include "family.h"
#include "matrix.h"
#include "random.h"

/* structure to hold estimation results */
typedef struct LME_struct {
  DIMS dm;          /* dimension data info */
  LENGTHS lengths;  /* lenghts and offsets object */
  FAMILY family;    /* family information */
  double
    *ZX,            /* model matrix */
    *y,             /* observed responses */
    *qraux,         /* auxiliar info for pre-decomposition */
    *settings,      /* settings for family object */
    *coef,          /* coefficients estimates */
    *theta,         /* scale parameters estimates */
    *scale,         /* scale estimate */
    *ranef,         /* random effects */
    *distances,     /* mahalanobis distances */
    *weights,       /* weights for heavy tailed distributions */
    *control,       /* control settings for estimation algorithm */
    *Delta,         /* relative precision factor */
    *Root;          /* triangular factors */
  int
    *dims,          /* dimensions passed from R */
    npar,           /* total length of the mixture parameters */
    maxIter,        /* maximun number of iterations */
    fixShape,       /* must estimate shape parameters? */
    ndraws,         /* independent draws for Monte Carlo integration */
    algorithm,      /* logical flag, EM = 0, nested EM = 1 */
    ncycles;        /* number of cycles for the nested EM algorithm */
  double
    tolerance;      /* convergence tolerance */
} LME_struct, *LME;

/* structure to hold fitted values */
typedef struct FITTED_struct {
  DIMS dm;          /* dimension data info */
  LENGTHS lengths;  /* lenghts and offsets object */
  double
    *ZX,            /* model matrix */
    *coef,          /* coefficient estimates */
    *ranef,         /* random effects */
    *conditional,   /* conditional fitted values */
    *marginal;      /* marginal fitted values */
} FITTED_struct, *FITTED;

/* structure to hold the covariance of coefficients */
typedef struct ACOV_struct {
  DIMS dm;          /* dimension data info */
  LENGTHS lengths;  /* lenghts and offsets object */
  FAMILY family;    /* family information */
  double
    *ZX,            /* model matrix */
    *Root,          /* triangular factors */
    *scale,         /* scale estimate */
    *control,       /* control settings */
    *acov;          /* coefficients covariance matrix */
  int
    ndraws;         /* independent draws for Monte Carlo integration */
} ACOV_struct, *ACOV;

/* estimation in lme under heavy tailed distributions (to be called by R) */
extern void lme_fit(double *, double *, double *, int *, int *, int *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *);
extern void lme_fitted(double *, int *, int *, int *, double *, double *, double *, double *);
extern void lme_acov(double *, int *, int *, int *, double *, double *, double *, double *, double *);

#endif /* HEAVY_LME_FIT_H */
