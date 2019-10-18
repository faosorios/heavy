/* ID: family.c, last updated 2019/08/02, F. Osorio */

#include "family.h"

/* static functions.. */
static double weight_normal();
static double weight_cauchy(double, double);
static double weight_student(double, double, double);
static double weight_slash(double, double, double);
static double weight_contaminated(double, double, double, double);

static double rand_weight_normal();
static double rand_weight_cauchy(double, double);
static double rand_weight_student(double, double, double);
static double rand_weight_slash(double, double, double);

static double negQfnc_student(double, void *);
static void update_df_student(DIMS, double *, double *, double *, double);
static double logwts_slash(double, double, double);
static void update_df_slash(DIMS, double *, double *, double *);

static double logLik_normal(DIMS, double *);
static double logLik_cauchy(DIMS, double *, double *);
static double logLik_student(DIMS, double *, double, double *);
static double logLik_slash(DIMS, double *, double, double *);
static double logLik_contaminated(DIMS, double *, double, double, double *);

static double acov_scale_normal();
static double acov_scale_cauchy(double);
static double acov_scale_student(double, double);
static double acov_scale_slash(double, double, int);
static double acov_scale_contaminated(double, double, double, int);
/* ..end declarations */

/* functions for dealing with 'family' objects */

FAMILY
family_init(double *settings)
{ /* constructor for a family object */
  FAMILY ans;

  ans = (FAMILY) Calloc(1, FAMILY_struct);
  ans->kind  = (int) settings[0];
  ans->npars = (int) settings[1];
  ans->nu = settings + 2;
  return ans;
}

void
family_free(FAMILY this)
{ /* destructor for a family object*/
  Free(this);
}

/* functions for computation of weights */

static double
weight_normal()
{ /* normal weight */
  return 1.;
}

static double
weight_cauchy(double length, double distance)
{ /* Cauchy weight */
  double wts;

  wts = (1. + length) / (1. + distance);
  return wts;
}

static double
weight_student(double length, double df, double distance)
{ /* Student-t weight */
  double wts;

  wts = (df + length) / (df + distance);
  return wts;
}

static double
weight_slash(double length, double df, double distance)
{ /* slash weight */
  int lower_tail = 1, log_p = 1;
  double shape, f, incr, wts;

  shape = df + .5 * length;
  if (shape > 178.) return 1.;
  f = (length + 2. * df) / distance;
  incr  = shape * log(distance / 2.) - .5 * distance - lgamma1p(shape);
  incr -= pgamma(1., shape, 2. / distance, lower_tail, log_p);
  wts = 1. - exp(incr);
  return (f * wts);
}

static double
weight_contaminated(double length, double epsilon, double vif, double distance)
{ /* contaminated normal weight */
  double f, wts;

  f = exp(.5 * (1. - vif) * distance);
  wts  = 1. - epsilon + epsilon * f * R_pow(vif, .5 * length + 1.);
  wts /= 1. - epsilon + epsilon * f * R_pow(vif, .5 * length);
  return wts;
}

double
do_weight(FAMILY family, double length, double distance)
{ /* weights dispatcher */
  double df, epsilon, vif, wts;

  switch (family->kind) {
    case NORMAL:
      wts = weight_normal();
      break;
    case CAUCHY:
      wts = weight_cauchy(length, distance);
      break;
    case STUDENT:
      df = (family->nu)[0];
      wts = weight_student(length, df, distance);
      break;
    case SLASH:
      df = (family->nu)[0];
      wts = weight_slash(length, df, distance);
      break;
    case CONTAMINATED:
      epsilon = (family->nu)[0];
      vif = (family->nu)[1];
      wts = weight_contaminated(length, epsilon, vif, distance);
      break;
    default:
      wts = weight_normal();
      break;
  }
  return wts;
}

static double
rand_weight_normal()
{ /* Normal: 'generation'? from the weights distribution (NOT to be used) */
  return 1.;
}

static double
rand_weight_cauchy(double length, double distance)
{ /* Cauchy: generation from the weights distribution */
  double val;

  val = rgamma(.5 * (1. + length), 2. / (1. + distance));
  return val;
}

static double
rand_weight_student(double length, double df, double distance)
{ /* Student-t: generation from the weights distribution */
  double val;

  val = rgamma(.5 * (df + length), 2. / (df + distance));
  return val;
}

static double
rand_weight_slash(double length, double df, double distance)
{ /* slash: generation from the weights distribution */
  double val;

  val = rng_tgamma_standard(df + .5 * length, .5 * distance);
  return val;
}

double
rand_weight(FAMILY family, double length, double distance)
{ /* weights dispatcher */
  double df, wts;

  switch (family->kind) {
    case NORMAL:
      wts = rand_weight_normal();
      break;
    case CAUCHY:
      wts = rand_weight_cauchy(length, distance);
      break;
    case STUDENT:
      df = (family->nu)[0];
      wts = rand_weight_student(length, df, distance);
      break;
    case SLASH:
      df = (family->nu)[0];
      wts = rand_weight_slash(length, df, distance);
      break;
    case CONTAMINATED:
      wts = rand_weight_normal(); /* FIXME */
      break;
    default:
      wts = rand_weight_normal();
      break;
  }
  return wts;
}

/* functions for the estimation of the mixture 'shape' parameters */

double
negQfnc_student(double df, void *pars)
{ /* for brent procedure */
  QTpars st = (QTpars) pars;
  DIMS dm = st->dm;
  double accum = 0.0, val;

  for (int i = 0; i < dm->n; i++) {
    accum += log((st->weights)[i]) - (st->weights)[i];
    accum += digamma(.5 * (st->df + (st->lengths)[i])) - log(.5 * (st->df + (st->lengths)[i]));
  }
  accum /= dm->n;

  /* compute Q-function for Student-t */
  val  = .5 * df * log(.5 * df) - lgammafn(.5 * df);
  val += .5 * df * accum;
  val *= dm->n;
  st->Qfnc = val;

  return -val;
}

static void
update_df_student(DIMS dm, double *df, double *weights, double *lengths, double tol)
{
  double conv, upper_df;
  const double c = (1. + sqrt(5.)) * .5;
  QTpars pars;

  pars = (QTpars) Calloc(1, QT_pars);

  /* constructs a Q-function object */
  pars->dm = dm;
  pars->weights = weights;
  pars->lengths = lengths;
  pars->df = *df;

  /* call optimizer */
  upper_df = *df;
  do {
    *df = brent(0., upper_df, negQfnc_student, pars, tol);
    conv = fabs(*df - upper_df);
    upper_df *= c;
  } while (conv < tol);

  Free(pars);
}

static double
logwts_slash(double length, double df, double distance)
{ /* slash log-weight */
  int lower_tail = 1, log_p = 0;
  double shape, incr, wts;

  shape = df + .5 * length;
  wts   = digamma(shape) - log(distance / 2.);
  incr  = pgamma_1st_derivative(1., shape, distance / 2.);
  incr /= pgamma(1., shape, 2. / distance, lower_tail, log_p);
  wts  += incr;
  return wts;
}

static void
update_df_slash(DIMS dm, double *df, double *distances, double *lengths)
{ /* update df for slash distribution */
  double accum = 0.0;

  for (int i = 0; i < dm->n; i++)
    accum += logwts_slash(lengths[i], *df, distances[i]);

  *df  = dm->n;
  *df /= -accum;
}

void
update_mixture(FAMILY family, DIMS dm, double *distances, double *lengths,
  double *weights, double tol)
{ /* update dispatcher */
  switch (family->kind) {
    case NORMAL:
      break;
    case CAUCHY:
      break;
    case STUDENT:
      update_df_student(dm, family->nu, weights, lengths, tol);
      break;
    case SLASH:
      update_df_slash(dm, family->nu, distances, lengths);
      break;
    case CONTAMINATED:
      break; /* FIXME */
    default:
      break;
  }
}

/*  functions for evaluation of the log-likelihood */

static double
logLik_normal(DIMS dm, double *distances)
{ /* gaussian log-likelihood */
  double accum = 0.0;

  for (int i = 0; i < dm->n; i++)
    accum += *distances++;
  return (-.5 * accum - dm->N * M_LN_SQRT_2PI);
}

static double
logLik_cauchy(DIMS dm, double *lengths, double *distances)
{ /* Cauchy log-likelihood */
  double accum = 0.0, p;

  for (int i = 0; i < dm->n; i++) {
    p = *lengths++;
    accum += lgammafn(.5 * (p + 1.)) - .5 * (p + 1.) * log1p(*distances++);
  }
  return (accum - (dm->N + dm->n) * M_LN_SQRT_PI);
}

static double
logLik_student(DIMS dm, double *lengths, double df, double *distances)
{ /* Student-t log-likelihood */
  double accum = 0.0, c, p;

  c = dm->n * lgammafn(.5 * df) + .5 * dm->N * (log(df) + 2. * M_LN_SQRT_PI);
  for (int i = 0; i < dm->n; i++) {
    p = *lengths++;
    accum += lgammafn(.5 * (df + p)) - .5 * (df + p) * log1p(*distances++ / df);
  }
  return (accum - c);
}

static double
logLik_slash(DIMS dm, double *lengths, double df, double *distances)
{ /* Slash log-likelihood */
  int lower_tail = 1, log_p = 1;
  double accum = 0.0, p, u, shape;

  for (int i = 0; i < dm->n; i++) {
    p = *lengths++;
    u = *distances++;
    shape = df + p / 2.;
    accum += lgammafn(shape) + shape * log(2. / u);
    accum += pgamma(1., shape, 2. / u, lower_tail, log_p);
  }
  return (accum + dm->n * log(df) - dm->N * M_LN_SQRT_2PI);
}

static double
logLik_contaminated(DIMS dm, double *lengths, double eps, double vif, double *distances)
{ /* contaminated-normal log-likelihood */
  double accum = 0.0, p, f, u;

  for(int i = 0; i < dm->n; i++) {
    p = *lengths++;
    u = *distances++;
    f = eps * pow(vif, .5 * p) * exp(-.5 * vif * u) + (1. - eps) * exp(-.5 * u);
    accum += log(f);
  }
  return (accum - dm->N * M_LN_SQRT_2PI);
}

double
logLik_kernel(FAMILY family, DIMS dm, double *lengths, double *distances)
{ /* logLik dispatcher */
  double df, eps, vif, ans;

  switch (family->kind) {
    case NORMAL:
      ans = logLik_normal(dm, distances);
      break;
    case CAUCHY:
      ans = logLik_cauchy(dm, lengths, distances);
      break;
    case STUDENT:
      df = (family->nu)[0];
      ans = logLik_student(dm, lengths, df, distances);
      break;
    case SLASH:
      df = (family->nu)[0];
      ans = logLik_slash(dm, lengths, df, distances);
      break;
    case CONTAMINATED:
      eps = (family->nu)[0];
      vif = (family->nu)[1];
      ans = logLik_contaminated(dm, lengths, eps, vif, distances);
      break;
    default:
      ans = logLik_normal(dm, distances);
      break;
  }
  return ans;
}

/* scale factor required for the Fisher information matrix */

static double
acov_scale_normal()
{ /* normal scale */
  return 1.;
}

static double
acov_scale_cauchy(double length)
{ /* Cauchy scale */
  double acov;

  acov = (length + 1.) / (length + 3.);
  return acov;
}

static double
acov_scale_student(double length, double df)
{ /* Student-t scale */
  double acov;

  acov = (df + length) / (df + length + 2.);
  return acov;
}

static double
acov_scale_slash(double length, double df, int ndraws)
{ /* slash scale */
  double accum = 0., u, w, *z;

  if (df > 30.)
    return 1.;
  z = (double *) Calloc(length, double);

  GetRNGstate();
  for (int i = 0; i < ndraws; i++) {
    rand_spherical_slash(z, df, 1, length);
    u = norm_sqr(z, 1, length);
    w = weight_slash(length, df, u);
    accum += SQR(w) * u;
  }
  PutRNGstate();
  accum /= ndraws;

  Free(z);
  return (accum / length);
}

static double
acov_scale_contaminated(double length, double epsilon, double vif, int ndraws)
{ /* contaminated normal scale */
  double accum = 0., u, w, *z;

  z = (double *) Calloc(length, double);

  GetRNGstate();
  for (int i = 0; i < ndraws; i++) {
    rand_spherical_contaminated(z, epsilon, vif, 1, length);
    u = norm_sqr(z, 1, length);
    w = weight_contaminated(length, epsilon, vif, u);
    accum += SQR(w) * u;
  }
  PutRNGstate();
  accum /= ndraws;

  Free(z);
  return (accum / length);
}

double
acov_scale(FAMILY family, double length, int ndraws)
{ /* scale factor for the Fisher information matrix */
  double df, epsilon, vif, ans;

  switch (family->kind) {
    case NORMAL:
      ans = acov_scale_normal();
      break;
    case CAUCHY:
      ans = acov_scale_cauchy(length);
      break;
    case STUDENT:
      df = (family->nu)[0];
      ans = acov_scale_student(length, df);
      break;
    case SLASH:
      df = (family->nu)[0];
      ans = acov_scale_slash(length, df, ndraws);
      break;
    case CONTAMINATED:
      epsilon = (family->nu)[0];
      vif = (family->nu)[1];
      ans = acov_scale_contaminated(length, epsilon, vif, ndraws);
      break;
    default:
      ans = acov_scale_normal();
      break;
  }
  return ans;
}
