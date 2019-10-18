/* ID: ps_fit.c, last updated 2019/08/05, F. Osorio */

#include "ps_fit.h"

/* static functions.. */
static DIMS dims(int *);
static void dims_free(DIMS);

static SPLINE ps_init(double *, double *, double *, int *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *);
static void ps_free(SPLINE);
static int ps_iterate(SPLINE);
static void ps_estimate(double *, double *, double *, DIMS, double *, double *, double *, double *, double *, double *, double *, double *, double *);
static int combined_iterate(SPLINE);
static void Mstep_and_WGCV(double *, double *, double *, DIMS, double *, double *, double *, double *, double *, double *, double *, double *, double *, double);
static double log_WGCV(double, void *);
static double plogLik(FAMILY, DIMS, double *, double *, double *, double *);
/* static functions.. */

static DIMS
dims(int *pdims)
{ /* dims object */
  DIMS dm;

  dm = (DIMS) Calloc(1, DIMS_struct);

  dm->N = (int) pdims[0];
  dm->n = dm->N;
  dm->p = (int) pdims[1];
  dm->deg = (int) pdims[2];
  dm->ord = (int) pdims[3];

  return dm;
}

static void
dims_free(DIMS this)
{ /* destructor for a dims object */
  Free(this);
}

static SPLINE
ps_init(double *y, double *b, double *half, int *pdims, double *settings, double *coef,
  double *scale, double *lambda, double *edf, double *gcv, double *pen, double *fitted,
  double *residuals, double *distances, double *weights, double *control)
{ /* constructor for a model object */
  SPLINE model;

  model = (SPLINE) Calloc(1, SPLINE_struct);
  model->dm = dims(pdims);
  model->family = family_init(settings);
  model->pdims = pdims;
  model->settings = settings;
  model->y = y;
  model->b = b;
  model->half = half;
  model->coef = coef;
  model->scale = scale;
  model->lambda = lambda;
  model->edf = edf;
  model->GCV = gcv;
  model->pen = pen;
  model->fitted = fitted;
  model->resid = residuals;
  model->distances = distances;
  model->weights = weights;
  model->control = control;
  model->maxIter = (int) control[0];
  model->tolerance = control[1];
  model->fixShape = (int) control[2];
  model->ndraws = (int) control[3];
  model->ncycles = (int) control[4];
  return model;
}

static void
ps_free(SPLINE this)
{ /* destructor for a model object */
  dims_free(this->dm);
  family_free(this->family);
  Free(this);
}

void
ps_fit(double *y, double *b, double *half, int *pdims, double *settings, double *coef,
  double *scale, double *lambda, double *edf, double *gcv, double *pen, double *fitted,
  double *residuals, double *distances, double *weights, double *logLik, double *control)
{ /* P-splines estimation */
  SPLINE model;

  model = ps_init(y, b, half, pdims, settings, coef, scale, lambda, edf, gcv, pen,
                  fitted, residuals, distances, weights, control);
  control[5] = ps_iterate(model);
  *logLik = plogLik(model->family, model->dm, model->distances, model->scale, model->lambda, model->pen);
  ps_free(model);
}

static int
ps_iterate(SPLINE model)
{ /* iteratively weighted P-splines estimation */
  DIMS dm = model->dm;
  int i, iter = 0;
  double conv, *lengths, logLik, newlogLik, *scale;

  /* initialization */
  lengths = (double *) Calloc(dm->n, double);
  for (i = 0; i < dm->n; i++)
    lengths[i] = 1.0;
  logLik = plogLik(model->family, model->dm, model->distances, model->scale, model->lambda, model->pen);

  /* main loop */
  repeat {
    /* E-step */
    scale = model->scale;
    for (i = 0; i < dm->n; i++) {
      (model->distances)[i] = SQR((model->resid)[i]) / *scale;
      (model->weights)[i] = do_weight(model->family, 1.0, (model->distances)[i]);
    }

    /* M-step */
    ps_estimate(model->y, model->b, model->half, model->dm, model->weights, model->coef,
                model->scale, model->lambda, model->GCV, model->edf, model->pen,
                model->fitted, model->resid);
    if (!model->fixShape)
      update_mixture(model->family, model->dm, model->distances, lengths, model->weights, model->tolerance);
    newlogLik = plogLik(model->family, model->dm, model->distances, model->scale, model->lambda, model->pen);

    iter++;

    /* eval convergence */
    conv = fabs((newlogLik - logLik) / (newlogLik + ABSTOL));
    if (conv < model->tolerance)
      break; /* successful completion */
    if (iter >= model->maxIter)
      break; /* maximum number of iterations exceeded */
    logLik = newlogLik;
  }
  Free(lengths);
  return iter;
}

static void
ps_estimate(double *y, double *b, double *half, DIMS dm, double *weights, double *coef,
  double *scale, double *lambda, double *gcv, double *edf, double *pen, double *fitted,
  double *resid)
{ /* weighted P-splines estimation */
  int i, j, qrows = dm->p - dm->ord, job, info = 0;
  double *a, *dummy = NULL, *u, *d, *v, *q, *r, *s, *rhs, *z, wts;
  double df = 0.0, div, GCV, PEN = 0.0, RSS, term;

  a    = (double *) Calloc(dm->p, double);
  u    = (double *) Calloc(dm->n * dm->p, double);
  d    = (double *) Calloc(dm->p, double);
  v    = (double *) Calloc(dm->p * dm->p, double);
  q    = (double *) Calloc(qrows * dm->p, double);
  r    = (double *) Calloc(dm->p, double);
  s    = (double *) Calloc(dm->p * dm->p, double);
  rhs  = (double *) Calloc(dm->p, double);
  z    = (double *) Calloc(dm->n, double);

  /* transform the model matrix and response */
  for (i = 0; i < dm->n; i++) {
    wts = sqrt(weights[i]);
    z[i] = wts * y[i];
    for (j = 0; j < dm->p; j++)
      u[i + j * dm->n] = wts * b[i + j * dm->n];
  }

  /* SVD of the 'model' matrix */
  job = 21; /* left singular vectors overwrite u */
  svd_decomp(u, dm->n, dm->n, dm->p, u, dm->n, d, v, dm->p, job, &info);
  if (info)
    error("1st call to svd_decomp in ps_estimate gave code %d", info);

  /* compute the right-hand side */
  crossprod(rhs, u, dm->n, dm->n, dm->p, z, dm->n, dm->n, 1);

  /* transform the half-penalty matrix */
  mult_mat(q, half, qrows, qrows, dm->p, v, dm->p, dm->p, dm->p);
  for (i = 0; i < qrows; i++) {
    for (j = 0; j < dm->p; j++)
      q[i + j * qrows] = q[i + j * qrows] / d[j];
  }

  /* SVD of the half-penalty matrix */
  job = 01; /* only compute singular values and right vectors */
  svd_decomp(q, qrows, qrows, dm->p, dummy, 0, r, s, dm->p, job, &info);
  if (info)
    error("2nd call to svd_decomp in ps_estimate gave code %d", info);

  /* update the right-hand side */
  crossprod(rhs, s, dm->p, dm->p, dm->p, rhs, dm->p, dm->p, 1);

  /* compute the coefficients and degrees of freedom */
  for (j = 0; j < dm->p; j++) {
    div  = 1.0 + *lambda * SQR(r[j]);
    df  += 1.0 / div;
    term = rhs[j] * r[j] / div;
    PEN += SQR(term);
    a[j] = rhs[j] / div;
  }
  mult_mat(a, s, dm->p, dm->p, dm->p, a, dm->p, dm->p, 1);

  /* compute weighted fitted values and residuals */
  mult_mat(fitted, u, dm->n, dm->n, dm->p, a, dm->p, dm->p, 1);
  for (i = 0; i < dm->n; i++)
    resid[i] = z[i] - fitted[i];

  /* compute GCV criterion */
  RSS  = norm_sqr(resid, 1, dm->n);
  GCV  = RSS / dm->n;
  GCV /= SQR(1.0 - df / dm->n);

  /* compute coefficients, unweighted fitted values and residuals */
  for (j = 0; j < dm->p; j++)
    coef[j] = a[j] / d[j];
  mult_mat(coef, v, dm->p, dm->p, dm->p, coef, dm->p, dm->p, 1);
  for (i = 0; i < dm->n; i++) {
    wts = sqrt(weights[i]);
    resid[i] /= wts;
    fitted[i] /= wts;
  }

  *edf = df;
  *gcv = log(GCV);
  *pen = PEN;
  *scale = (RSS + *lambda * PEN) / dm->n;

  Free(a); Free(u); Free(d); Free(v); Free(q); Free(r); Free(s); Free(rhs); Free(z);
}

void
ps_combined(double *y, double *b, double *half, int *pdims, double *settings, double *coef,
  double *scale, double *lambda, double *edf, double *gcv, double *pen, double *fitted,
  double *residuals, double *distances, double *weights, double *logLik, double *control)
{ /* P-splines estimation with selection of the smoothing parameter */
  SPLINE model;

  model = ps_init(y, b, half, pdims, settings, coef, scale, lambda, edf, gcv, pen, fitted,
                  residuals, distances, weights, control);
  control[5] = combined_iterate(model);
  *logLik = plogLik(model->family, model->dm, model->distances, model->scale, model->lambda, model->pen);
  ps_free(model);
}

static int
combined_iterate(SPLINE model)
{ /* iteratively weighted P-splines estimation */
  DIMS dm = model->dm;
  int i, iter = 0;
  double conv, *lengths, logLik, newlogLik, *scale;

  /* initialization */
  lengths = (double *) Calloc(dm->n, double);
  for (i = 0; i < dm->n; i++)
    lengths[i] = 1.0;
  logLik = plogLik(model->family, model->dm, model->distances, model->scale, model->lambda, model->pen);

  /* main loop */
  repeat {
    /* E-step */
    scale = model->scale;
    for (i = 0; i < dm->n; i++) {
      (model->distances)[i] = SQR((model->resid)[i]) / *scale;
      (model->weights)[i] = do_weight(model->family, 1.0, (model->distances)[i]);
    }

    /* combined M-step and WGCV search */
    Mstep_and_WGCV(model->y, model->b, model->half, model->dm, model->weights, model->coef,
                   model->scale, model->lambda, model->GCV, model->edf, model->pen,
                   model->fitted, model->resid, model->tolerance);
    if (!model->fixShape)
      update_mixture(model->family, model->dm, model->distances, lengths, model->weights, model->tolerance);
    newlogLik = plogLik(model->family, model->dm, model->distances, model->scale, model->lambda, model->pen);

    iter++;

    /* eval convergence */
    conv = fabs((newlogLik - logLik) / (newlogLik + ABSTOL));
    if (conv < model->tolerance)
      break; /* successful completion */
    if (iter >= model->maxIter)
      break; /* maximum number of iterations exceeded */
    logLik = newlogLik;
  }
  Free(lengths);
  return iter;
}

static void
Mstep_and_WGCV(double *y, double *b, double *half, DIMS dm, double *weights, double *coef,
  double *scale, double *lambda, double *gcv, double *edf, double *pen, double *fitted,
  double *resid, double tol)
{ /* both weighted B-splines and GCV iterations */
  int i, j, qrows = dm->p - dm->ord, job, info = 0;
  double *a, *dummy = NULL, *u, *d, *v, *q, *r, *s, *rhs, *z;
  double conv, upper_lambda, wts;
  const double c = (1.0 + sqrt(5.0)) * 0.5;
  GCVpars pars;

  a    = (double *) Calloc(dm->p, double);
  u    = (double *) Calloc(dm->n * dm->p, double);
  d    = (double *) Calloc(dm->p, double);
  v    = (double *) Calloc(dm->p * dm->p, double);
  q    = (double *) Calloc(qrows * dm->p, double);
  r    = (double *) Calloc(dm->p, double);
  s    = (double *) Calloc(dm->p * dm->p, double);
  rhs  = (double *) Calloc(dm->p, double);
  z    = (double *) Calloc(dm->n, double);
  pars = (GCVpars) Calloc(1, GCV_pars);

  /* transform the model matrix and response */
  for (i = 0; i < dm->n; i++) {
    wts = sqrt(weights[i]);
    z[i] = wts * y[i];
    for (j = 0; j < dm->p; j++)
      u[i + j * dm->n] = wts * b[i + j * dm->n];
  }

  /* SVD of the 'model' matrix */
  job = 21; /* left singular vectors overwrite u */
  svd_decomp(u, dm->n, dm->n, dm->p, u, dm->n, d, v, dm->p, job, &info);
  if (info)
    error("1st call to svd_decomp in Mstep_and_WGCV gave code %d", info);

  /* compute the right-hand side */
  crossprod(rhs, u, dm->n, dm->n, dm->p, z, dm->n, dm->n, 1);

  /* transform the half-penalty matrix */
  mult_mat(q, half, qrows, qrows, dm->p, v, dm->p, dm->p, dm->p);
  for (i = 0; i < qrows; i++) {
    for (j = 0; j < dm->p; j++)
      q[i + j * qrows] = q[i + j * qrows] / d[j];
  }

  /* SVD of the half-penalty matrix */
  job = 01; /* only compute singular values and right vectors */
  svd_decomp(q, qrows, qrows, dm->p, dummy, 0, r, s, dm->p, job, &info);
  if (info)
    error("2nd call to svd_decomp in Mstep_and_WGCV gave code %d", info);

  /* update the right-hand side */
  crossprod(rhs, s, dm->p, dm->p, dm->p, rhs, dm->p, dm->p, 1);

  /* constructs a GCV object */
  pars->dm = dm;
  pars->u = u;
  pars->r = r;
  pars->s = s;
  pars->z = z;
  pars->rhs = rhs;
  pars->a = a;
  pars->fitted = fitted;
  pars->resid  = resid;

  /* call optimizer */
  upper_lambda = *lambda;
  do {
    *lambda = brent(0., upper_lambda, log_WGCV, pars, tol);
    conv = fabs(*lambda - upper_lambda);
    upper_lambda *= c;
  } while (conv < tol);

  /* compute coefficients, unweighted fitted values and residuals */
  for (j = 0; j < dm->p; j++)
    coef[j] = a[j] / d[j];
  mult_mat(coef, v, dm->p, dm->p, dm->p, coef, dm->p, dm->p, 1);
  for (i = 0; i < dm->n; i++) {
    wts = sqrt(weights[i]);
    resid[i] /= wts;
    fitted[i] /= wts;
  }

  *edf = pars->edf;
  *gcv = pars->WGCV;
  *pen = pars->pen;
  *scale = (pars->RSS + *lambda * pars->pen) / dm->n;

  Free(a); Free(u); Free(d); Free(v); Free(q); Free(r); Free(s); Free(rhs); Free(z); Free(pars);
}

double
log_WGCV(double lambda, void *pars)
{ /* for brent's procedure */
  GCVpars st = (GCVpars) pars;
  DIMS dm = st->dm;
  double div, edf = 0.0, PEN = 0.0, s2, term, val;

  /* compute the coefficients and degrees of freedom */
  for (int j = 0; j < dm->p; j++) {
    div  = 1.0 + lambda * SQR((st->r)[j]);
    edf += 1.0 / div;
    term = (st->rhs)[j] * (st->r)[j] / div;
    PEN += SQR(term);
    (st->a)[j] = (st->rhs)[j] / div;
  }
  mult_mat(st->a, st->s, dm->p, dm->p, dm->p, st->a, dm->p, dm->p, 1);

  /* compute weighted fitted values and residuals */
  mult_mat(st->fitted, st->u, dm->n, dm->n, dm->p, st->a, dm->p, dm->p, 1);
  for (int i = 0; i < dm->n; i++)
    (st->resid)[i] = (st->z)[i] - (st->fitted)[i];

  /* compute log-WGCV criteria */
  st->RSS = norm_sqr(st->resid, 1, dm->n);
  s2 = st->RSS / (dm->n - edf);
  val = log(s2) - log(1.0 - edf / dm->n);
  st->edf = edf;
  st->pen = PEN;
  st->WGCV = val;

  return val;
}

static double
plogLik(FAMILY family, DIMS dm, double *distances, double *scale, double *lambda, double *pen)
{ /* evaluate the penalized log-likelihood function */
  double krnl, *lengths;

  lengths = (double *) Calloc(dm->n, double);
  for (int i = 0; i < dm->n; i++)
    lengths[i] = 1.0;
  krnl = logLik_kernel(family, dm, lengths, distances);
  return (-0.5 * dm->n * log(*scale) + krnl - 0.5 * *lambda * *pen / *scale);
}
