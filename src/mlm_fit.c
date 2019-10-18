/* ID: mlm_fit.c, last updated 2019/08/05, F. Osorio */

#include "lm_fit.h"

/* static functions.. */
static DIMS dims(int *);
static void dims_free(DIMS);

static LM mlm_init(double *, double *, int *, double *, double *, double *, double *, double *, double *, double *, double *, double *);
static void mlm_free(LM);
static int IRLS(double *, double *, DIMS, FAMILY, double *, double *, double *, double *, double *, double *, int, double, int);
static void IRLS_increment(double *, double *, DIMS, double *, double *, double *, double *, double *);
static double mahalanobis(double *, int, double *, int, double *, double *);

static double mlm_logLik(FAMILY, DIMS, double *, double *);
static void mlm_acov(FAMILY, DIMS, double *, int, double *);
/* ..end declarations */

void
mlm_fit(double *y, double *x, int *pdims, double *settings, double *coef, double *scale,
  double *fitted, double *resid, double *distances, double *weights, double *logLik,
  double *acov, double *control)
{ /* fitter for multivariate linear regression under heavy-tailed errors */
  LM model;

  model = mlm_init(y, x, pdims, settings, coef, scale, fitted, resid, distances, weights, acov, control);
  control[4] = (double) IRLS(model->y, model->x, model->dm, model->family, model->coef,
                             model->scale, model->resid, model->fitted, model->distances,
                             model->weights, model->maxIter, model->tolerance, model->fixShape);
  *logLik = mlm_logLik(model->family, model->dm, model->distances, model->scale);
  mlm_acov(model->family, model->dm, model->x, model->ndraws, model->acov);
  mlm_free(model);
}

static DIMS
dims(int *pdims)
{ /* dims object for multivariate linear models */
  DIMS ans;

  ans = (DIMS) Calloc(1, DIMS_struct);
  ans->N  = (int) pdims[0];
  ans->n  = ans->N;
  ans->p  = (int) pdims[1];
  ans->ny = (int) pdims[2];
  return ans;
}

static void
dims_free(DIMS this)
{ /* destructor for a dims object */
  Free(this);
}

static LM
mlm_init(double *y, double *x, int *pdims, double *settings, double *coef,
  double *scale, double *fitted, double *resid, double *distances, double *weights,
  double *acov, double *control)
{ /* constructor for a multivariate linear model object */
  LM model;

  model = (LM) Calloc(1, LM_struct);
  model->dm = dims(pdims);
  model->settings = settings;
  model->family = family_init(settings);
  model->y = y;
  model->x = x;
  model->coef = coef;
  model->scale = scale;
  model->fitted = fitted;
  model->resid = resid;
  model->distances = distances;
  model->weights = weights;
  model->acov = acov;
  model->control = control;
  model->maxIter = (int) control[0];
  model->tolerance = control[1];
  model->fixShape = (int) control[2];
  model->ndraws = (int) control[3];
  return model;
}

static void
mlm_free(LM this)
{ /* destructor for a model object */
  dims_free(this->dm);
  family_free(this->family);
  Free(this);
}

static int
IRLS(double *y, double *x, DIMS dm, FAMILY family, double *coef, double *scale,
  double *resid, double *fitted, double *distances, double *weights, int maxit,
  double tolerance, int fixShape)
{ /* iteratively reweighted LS algorithm */
  int iter = 0, rdf = dm->n - dm->p;
  double conv, tol = R_pow(tolerance, 2.0/3.0), RSS, newRSS;
  double *lengths, *u, *z, *working;

  /* initialization */
  lengths = (double *) Calloc(dm->n, double);
  u       = (double *) Calloc(dm->ny, double);
  z       = (double *) Calloc(dm->p, double);
  working = (double *) Calloc(dm->n * dm->ny, double);
  for (int i = 0; i < dm->n; i++)
    lengths[i] = (double) dm->ny;
  RSS = (double) dm->n * dm->ny;

  /* main loop */
  repeat {
    /* E-step */
    for (int i = 0; i < dm->n; i++) {
      copy_vec(u, 1, y + i, dm->n, dm->ny);
      copy_vec(z, 1, x + i, dm->n, dm->p);
      distances[i] = mahalanobis(u, dm->ny, z, dm->p, coef, scale);
      weights[i] = do_weight(family, lengths[i], distances[i]);
    }

    /* M-step */
    IRLS_increment(y, x, dm, coef, resid, fitted, weights, working);
    crossprod(scale, working + dm->p, dm->n, rdf, dm->ny, working + dm->p, dm->n, rdf, dm->ny);
    scale_mat(scale, dm->ny, (double) 1.0/dm->n, scale, dm->ny, dm->ny, dm->ny);
    if (!fixShape)
      update_mixture(family, dm, distances, lengths, weights, tol);

    iter++;

    /* eval convergence */
    newRSS = dot_product(weights, 1, distances, 1, dm->n);
    conv = fabs((newRSS - RSS) / (newRSS + ABSTOL));
    if (conv < tolerance)
      break; /* successful completion */
    if (iter >= maxit)
      break; /* maximum number of iterations exceeded */
    RSS = newRSS;
  }
  Free(lengths); Free(u); Free(z); Free(working);
  return iter;
}

static void
IRLS_increment(double *y, double *x, DIMS dm, double *coef, double *resid, double *fitted,
  double *weights, double *working)
{ /* increment for direction search in IRLS */
  int i, j, info = 0;
  double one = 1.0, stepsize = 1.0, wts, *incr, *qraux, *u, *z;
  char *side = "L", *uplo = "U", *diag = "N", *notrans = "N";
  QRStruct qr;

  /* initialization */
  incr  = (double *) Calloc(dm->p * dm->ny, double);
  qraux = (double *) Calloc(dm->p, double);
  u     = (double *) Calloc(dm->n * dm->ny, double);
  z     = (double *) Calloc(dm->n * dm->p, double);

  /* transformed model matrix and working residuals */
  for (i = 0; i < dm->n; i++) {
    wts = sqrt(weights[i]);
    for (j = 0; j < dm->p; j++)
      z[i + j * dm->n] = wts * x[i + j * dm->n];
    for (j = 0; j < dm->ny; j++)
      working[i + j * dm->n] = wts * resid[i + j * dm->n];
  }

  /* solve the transformed LS-problem */
  qr = QR_decomp(z, dm->n, dm->n, dm->p, qraux);
  QR_coef(qr, incr, u, dm->n, dm->ny, &info);
  if (info)
    error("QR_coef in IRLS_increment gave code %d", info);

  /* update coefficients */
  add_mat(coef, dm->p, stepsize, incr, dm->p, dm->p, dm->ny);

  /* compute fitted values */
  setzero(u, dm->n, dm->n, dm->ny);
  copy_mat(u, dm->n, coef, dm->p, dm->p, dm->ny);
  F77_CALL(dtrmm)(side, uplo, notrans, diag, &(dm->p), &(dm->ny), &one, z, &(dm->n), u, &(dm->n));
  QR_qy(qr, fitted, u, dm->n, dm->ny, &info);
  if (info)
    error("qr_qy in IRLS_increment gave code %d", info);

  /* un-weighted fitted values and residuals */
  for (i = 0; i < dm->n; i++) {
    wts = sqrt(weights[i]);
    for (j = 0; j < dm->ny; j++) {
      fitted[i + j * dm->n] /= wts;
      resid[i + j * dm->n] = y[i + j * dm->n] - fitted[i + j * dm->n];
    }
  }

  QR_free(qr); Free(incr); Free(qraux); Free(u); Free(z);
}

static double
mahalanobis(double *y, int ny, double *x, int p, double *coef, double *Scatter)
{ /* Mahalanobis distances */
  double ans, *center, *z, *Root;
  int info = 0, job = 0;

  center = (double *) Calloc(ny, double);
  z      = (double *) Calloc(ny, double);
  Root   = (double *) Calloc(ny * ny, double);

  copy_mat(Root, ny, Scatter, ny, ny, ny);
  chol_decomp(Root, ny, ny, job, &info);
  if (info)
    error("chol_decomp in mahalanobis gave code %d", info);

  Memcpy(z, y, ny);
  crossprod(center, coef, p, p, ny, x, p, p, 1);
  ax_plus_y(-1.0, center, 1, z, 1, ny);
  backsolve(Root, ny, ny, z, ny, 1, job, &info);
  if (info)
    error("backsolve in mahalanobis gave code %d", info);
  ans = norm_sqr(z, 1, ny);

  Free(center); Free(z); Free(Root);
  return ans;
}

static double
mlm_logLik(FAMILY family, DIMS dm, double *distances, double *Scatter)
{ /* evaluate the log-likelihood function for multivariate linear models */
  double ans = 0.0, *lengths, *Root;
  int info = 0, job = 0;

  lengths = (double *) Calloc(dm->n, double);
  for (int i = 0; i < dm->n; i++)
    lengths[i] = (double) dm->ny;
  Root = (double *) Calloc(SQR(dm->ny), double);
  copy_mat(Root, dm->ny, Scatter, dm->ny, dm->ny, dm->ny);
  chol_decomp(Root, dm->ny, dm->ny, job, &info);
  if (info)
    error("chol_decomp in mlm_logLik gave code %d", info);
  ans -= dm->n * logAbsDet(Root, dm->ny, dm->ny);
  ans += logLik_kernel(family, dm, lengths, distances);
  Free(lengths); Free(Root);
  return ans;
}

static void
mlm_acov(FAMILY family, DIMS dm, double *x, int ndraws, double *acov)
{ /* evaluate the Fisher information matrix for the coefficients */
  int info = 0, job = 1;
  double factor, *qraux, *R;
  QRStruct qr;

  /* initialization */
  qraux = (double *) Calloc(dm->p, double);
  R     = (double *) Calloc(SQR(dm->p), double);

  /* unscaled Fisher information matrix */
  qr = QR_decomp(x, dm->n, dm->n, dm->p, qraux);
  QR_store_R(qr, R, dm->p);
  invert_triangular(R, dm->p, dm->p, job, &info);
  if (info)
    error("invert_triangular in mlm_acov gave code %d", info);
  outerprod(acov, R, dm->p, dm->p, dm->p, R, dm->p, dm->p, dm->p);

  /* scaling */
  factor = 1. / acov_scale(family, (double) dm->ny, ndraws);
  scale_mat(acov, dm->p, factor, acov, dm->p, dm->p, dm->p);

  QR_free(qr); Free(qraux); Free(R);
}
