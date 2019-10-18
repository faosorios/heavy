/* ID: family.h, last updated 2019/08/02, F. Osorio */

#ifndef HEAVY_FAMILY_H
#define HEAVY_FAMILY_H

#include "base.h"
#include "distn.h"
#include "optim.h"
#include "random.h"
#include "specfun.h"

/* available families */
typedef enum {
  NORMAL,
  CAUCHY,
  STUDENT,
  SLASH,
  CONTAMINATED
} classes;

/* heavy tailed family structure */
typedef struct FAMILY_struct {
  classes kind;   /* family kind */
  int npars;      /* number of parameters in 'family' */
  double *nu;     /* parameter vector */
} FAMILY_struct, *FAMILY;

/* Q-function info required for the degrees of freedom estimation */
typedef struct QT_pars {
  DIMS dm;
  double df, Qfnc;
  double *lengths, *weights;
} QT_pars, *QTpars;

/* functions for dealing with 'family' objects */
extern FAMILY family_init(double *);
extern void family_free(FAMILY);

/* routines for computation of weights */
extern double do_weight(FAMILY, double, double);
extern double rand_weight(FAMILY, double, double);
extern void update_mixture(FAMILY, DIMS, double *, double *, double *, double);

/*  functions for evaluation of the log-likelihood */
extern double logLik_kernel(FAMILY, DIMS, double *, double *);

/* scale factor for the Fisher information matrix */
extern double acov_scale(FAMILY, double, int);

#endif /* HEAVY_FAMILY_H */
