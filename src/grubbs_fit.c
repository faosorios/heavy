/* ID: grubbs_fit.c, last updated 2019/08/23, F. Osorio */

#include "grubbs_fit.h"

/* static functions.. */
static DIMS dims(int *);
static void dims_free(DIMS);

static GRUBBS grubbs_init(double *, int *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *);
static void grubbs_free(GRUBBS);
static int grubbs_iterate(GRUBBS);
static void grubbs_Estep(double *, DIMS, FAMILY, double *, double *, double *, double *, double *, double *);
static void update_mu(double *, DIMS, double *, double *, double *, double *);
static void update_phi_and_scale(double *, DIMS, double *, double *, double *, double *, double *, double *);

static double grubbs_logLik(FAMILY, DIMS, double *, double *, double *);
static void grubbs_acov(FAMILY, DIMS, double *, double *, double *, int, double *);
/* ..end declarations */

void
grubbs_fit(double *y, int *pdims, double *settings, double *mu, double *phi, double *scale,
  double *z, double *distances, double *weights, double *residuals, double *logLik, double *acov,
  double *control)
{ /* Grubbs' model estimation procedure under heavy-tailed errors */
  GRUBBS model;

  model = grubbs_init(y, pdims, settings, mu, phi, scale, z, distances, weights, residuals, acov, control);
  control[4] = grubbs_iterate(model);
  *logLik = grubbs_logLik(model->family, model->dm, model->distances, model->scale, model->phi);
  grubbs_acov(model->family, model->dm, model->mu, model->scale, model->phi, model->ndraws, model->acov);
  grubbs_free(model);
}

static DIMS
dims(int *pdims)
{ /* dims object for Grubbs' model */
  DIMS ans;

  ans = (DIMS) Calloc(1, DIMS_struct);
  ans->N = (int) pdims[0];
  ans->n = ans->N;
  ans->p = (int) pdims[1];
  return ans;
}

static void
dims_free(DIMS this)
{ /* destructor for a dims object */
  Free(this);
}

static GRUBBS
grubbs_init(double *y, int *pdims, double *settings, double *mu, double *phi, double *scale,
  double *z, double *distances, double *weights, double *residuals, double *acov, double *control)
{ /* constructor for a Grubbs object */
  GRUBBS model;

  model = (GRUBBS) Calloc(1, GRUBBS_struct);
  model->dm = dims(pdims);
  model->settings = settings;
  model->family = family_init(settings);
  model->y = y;
  model->mu = mu;
  model->phi = phi;
  model->scale = scale;
  model->z = z;
  model->distances = distances;
  model->weights = weights;
  model->residuals = residuals;
  model->acov = acov;
  model->control = control;
  model->maxIter = (int) control[0];
  model->tolerance = control[1];
  model->fixShape = (int) control[2];
  model->ndraws = (int) control[3];
  return model;
}

static void
grubbs_free(GRUBBS this)
{ /* destructor for a Grubbs object */
  dims_free(this->dm);
  family_free(this->family);
  Free(this);
}

static int
grubbs_iterate(GRUBBS model)
{
  int iter = 0;
  double conv, tol = R_pow(model->tolerance, 2.0/3.0), *lengths, RSS, newRSS;

  /* initialization */
  lengths = (double *) Calloc(model->dm->n, double);
  for (int i = 0; i < model->dm->n; i++)
    lengths[i] = (double) model->dm->p;
  RSS = (double) model->dm->n * model->dm->p;

  /* main loop */
  repeat {
    /* E-step */
    grubbs_Estep(model->y, model->dm, model->family, model->mu, model->phi, model->scale,
                 model->z, model->distances, model->weights);

    /* CM-steps */
    update_mu(model->y, model->dm, model->mu, model->phi, model->z, model->weights);
    update_phi_and_scale(model->y, model->dm, model->mu, model->phi, model->scale,
                         model->z, model->weights, model->residuals);
    if (!(model->fixShape))
      update_mixture(model->family, model->dm, model->distances, lengths, model->weights, tol);
    newRSS = dot_product(model->weights, 1, model->distances, 1, model->dm->n);

    iter++;

    /* eval convergence */
    conv = fabs((newRSS - RSS) / (newRSS + ABSTOL));
    if (conv < model->tolerance)
      break; /* successful completion */
    if (iter >= model->maxIter)
      break; /* maximum number of iterations exceeded */
    RSS = newRSS;
  }
  Free(lengths);
  return iter;
}

static void
grubbs_Estep(double *y, DIMS dm, FAMILY family, double *mu, double *phi, double *scale,
  double *z, double *distances, double *weights)

{
  double accum = 0.0, tau, zhat;

  /* computation of tau */
  tau = *scale;
  for (int j = 0; j < dm->p; j++)
    accum += 1.0 / phi[j];
  tau /= 1.0 + *scale * accum;

  /* computation of z_i estimates and distances */
  for (int i = 0; i < dm->n; i++) {
    double *rsp = y + i * dm->p;

    zhat = 0.0;
    for (int j = 0; j < dm->p; j++)
      zhat += (rsp[j] - mu[j]) / phi[j];
    zhat *= tau;

    accum = 0.0;
    for (int j = 0; j < dm->p; j++)
      accum += SQR(rsp[j] - mu[j] - zhat) / phi[j];

    z[i] = zhat;
    distances[i] = accum + SQR(zhat) / *scale;
    weights[i] = do_weight(family, (double) dm->p, distances[i]);
  }
}

static void
update_mu(double *y, DIMS dm, double *mu, double *phi, double *z, double *weights)
{
  double accum = 0.0, wts, zhat;

  /* initialization */
  setzero(mu, 0, dm->p, 1);

  /* compute mu estimate */
  for (int i = 0; i < dm->n; i++) {
    double *rsp = y + i * dm->p;
    wts  = weights[i];
    zhat = z[i];

    for (int j = 0; j < dm->p; j++)
      mu[j] += wts * (rsp[j] - zhat);
    accum += wts;
  }

  for (int j = 0; j < dm->p; j++)
    mu[j] /= accum;
}

static void
update_phi_and_scale(double *y, DIMS dm, double *mu, double *phi, double *scale,
  double *z, double *weights, double *residuals)
{
  double accum = 0.0, tau, wts, zhat;

  /* computation of tau */
  tau = *scale;
  for (int j = 0; j < dm->p; j++)
    accum += 1.0 / phi[j];
  tau /= 1.0 + *scale * accum;

  /* compute phi and scale estimates */
  setzero(phi, 0, dm->p, 1);
  *scale = 0.0;
  for (int i = 0; i < dm->n; i++) {
    double *rsp = y + i * dm->p;
    wts  = weights[i];
    zhat = z[i];

    for (int j = 0; j < dm->p; j++) {
      phi[j] += wts * SQR(rsp[j] - mu[j] - zhat);
      *scale += wts * SQR(zhat);
      residuals[j] = rsp[j] - mu[j];
    }
    residuals += dm->p;
  }

  *scale /= dm->n;
  *scale += tau;
  for (int j = 0; j < dm->p; j++)
    phi[j] = tau + phi[j] / dm->n;
}

static double
grubbs_logLik(FAMILY family, DIMS dm, double *distances, double *scale, double *phi)
{ /* evaluate the log-likelihood function for Grubbs' model */
  double accum = 0.0, krnl, logdet = 0.0, *lengths;

  /* initialization */
  lengths = (double *) Calloc(dm->n, double);
  for (int i = 0; i < dm->n; i++)
    lengths[i] = 1.0;

  /* computation of log(det) */
  for (int j = 0; j < dm->p; j++) {
    accum  += 1.0 / phi[j];
    logdet -= log(phi[j]);
  }
  logdet += log1p(*scale * accum);

  /* computation of log-likelihood */
  krnl = logLik_kernel(family, dm, lengths, distances);
  Free(lengths);
  return (krnl - .5 * dm->n * logdet);
}

static void
grubbs_acov(FAMILY family, DIMS dm, double *mu, double *scale, double *phi, int ndraws, double *acov)
{ /* evaluate the Fisher information matrix of center parameters */
  double factor, *ones;

  /* initialization */
  ones = (double *) Calloc(dm->p, double);
  for (int j = 0; j < dm->p; j++) {
    ones[j] = 1.0;
    acov[j * (dm->p + 1)] = phi[j];
  }

  rank1_update(acov, dm->p, dm->p, dm->p, *scale, ones, ones);
  factor = 0.25 * dm->p / acov_scale(family, dm->p, ndraws);
  scale_mat(acov, dm->p, factor, acov, dm->p, dm->p, dm->p);
  Free(ones);
}
