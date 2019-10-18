/* ID: lme_fit.c, last updated 2019/08/08, F. Osorio */

#include "lme_fit.h"

/* static functions.. */
static DIMS dims(int *);
static void dims_free(DIMS);

static LENGTHS setLengths(DIMS, int *, int *);
static void lengths_free(LENGTHS);

static LME lme_init(double *, double *, double *, int *, int *, int *, double *, double *, double *, double *, double *, double *, double *, double *, double *);
static void lme_free(LME);
static void pre_decomp(double *, double *, double *, DIMS, LENGTHS);

static int lme_iterate(LME);
static void internal_EMcycle(LME);
static void internal_Estep(LME);
static void outer_Estep(LME);

static void update_coef(LME);
static void update_scale(LME);
static void update_theta(LME);
static void update_nu(LME);
static void relative_precision(double *, int, double *, double *);
static void append_decomp(double *, double *, int, int, int, double *, int, double *, int);
static void random_effects(double *, int, DIMS, double *, double *);
static double mahalanobis(double *, int, DIMS, double *, double *);
static void working_response(double *, double *, int, DIMS, double *, double *);

static double lme_logLik(LME);

static FITTED lme_fitted_init(double *, int *, int *, int *, double *, double *, double *, double *);
static void lme_fitted_free(FITTED);
static void lme_fitted_values(FITTED);

static ACOV lme_acov_init(double *, int *, int *, int *, double *, double *, double *, double *, double *);
static void lme_acov_free(ACOV);
static void lme_acov_coef(ACOV);
/* ..end declarations */

void
lme_fit(double *ZX, double *y, double *qraux, int *pdims, int *lengths, int *DcLen,
  double *settings, double *coef, double *theta, double *scale, double *ranef,
  double *Root, double *distances, double *weights, double *logLik, double *control)
{ /* fitter for linear mixed-effects models under heavy-tailed errors */
  LME model;

  model = lme_init(ZX, y, qraux, pdims, lengths, DcLen, settings, coef, theta, scale,
                   ranef, Root, distances, weights, control);
  pre_decomp(model->ZX, model->y, model->qraux, model->dm, model->lengths);
  control[6] = (double) lme_iterate(model);
  *logLik = lme_logLik(model);
  lme_free(model);
}

static DIMS
dims(int *pdims)
{ /* dims object for LME models */
  DIMS dm;

  dm = (DIMS) Calloc(1, DIMS_struct);

  dm->n = (int) pdims[0];
  dm->q = (int) pdims[1];
  dm->p = (int) pdims[2];
  dm->N = (int) pdims[3];
  dm->ZXrows = (int) pdims[4];
  dm->ZXcols = (int) pdims[5];
  dm->DcRows = (int) pdims[6];
  return dm;
}

static void
dims_free(DIMS this)
{ /* destructor for a dims object */
  Free(this);
}

static LENGTHS
setLengths(DIMS dm, int *glen, int *DcLen)
{ /* constructor for a lenghts object */
  LENGTHS len;
  int accum1 = 0, accum2 = 0, accum3 = 0;

  len = (LENGTHS) Calloc(1, LENGTHS_struct);
  len->glen  = glen;
  len->DcLen = DcLen;
  len->offsets = (int *) Calloc(dm->n, int);
  len->ZXlen   = (int *) Calloc(dm->n, int);
  len->ZXoff   = (int *) Calloc(dm->n, int);
  len->DcOff   = (int *) Calloc(dm->n, int);
  (len->ZXlen)[0] = (len->glen)[0] * (dm->ZXcols + 1);
  for (int i = 1; i < dm->n; i++) {
    (len->ZXlen)[i] = (len->glen)[i] * (dm->ZXcols + 1);
    accum1 += (len->glen)[i-1];
    (len->offsets)[i] = accum1;
    accum2 += (len->ZXlen)[i-1];
    (len->ZXoff)[i] = accum2;
    accum3 += (len->DcLen)[i-1];
    (len->DcOff)[i] = accum3;
  }
  return len;
}

static void
lengths_free(LENGTHS this)
{ /* destructor for a lenghts object */
  Free(this->offsets);
  Free(this->ZXlen);
  Free(this->ZXoff);
  Free(this->DcOff);
  Free(this);
}

static LME
lme_init(double *ZX, double *y, double *qraux, int *pdims, int *lengths, int *DcLen,
  double *settings, double *coef, double *theta, double *scale, double *ranef,
  double *Root, double *distances, double *weights, double *control)
{ /* constructor for a lme object */
  LME model;
  int qsq;

  model = (LME) Calloc(1, LME_struct);
  model->ZX = ZX;
  model->y = y;
  model->qraux = qraux;
  model->dims = pdims;
  model->dm = dims(pdims);
  model->lengths = setLengths(model->dm, lengths, DcLen);
  model->settings = settings;
  model->family = family_init(settings);
  model->coef = coef;
  model->theta = theta;
  model->scale = scale;
  model->ranef = ranef;
  model->Root = Root;
  model->distances = distances;
  model->weights = weights;
  model->control = control;
  /* some definitions */
  qsq = SQR(model->dm->q);
  model->npar = model->dm->p + qsq + 1;
  model->Delta = (double *) Calloc(qsq, double);
  model->maxIter = (int) control[0];
  model->tolerance = control[1];
  model->fixShape = (int) control[2];
  model->ndraws = (int) control[3];
  model->algorithm = (int) control[4];
  model->ncycles = (int) control[5];
  return model;
}

static void
lme_free(LME this)
{ /* destructor for a lme object */
  dims_free(this->dm);
  lengths_free(this->lengths);
  family_free(this->family);
  Free(this->Delta);
  Free(this);
}

static void
pre_decomp(double *ZX, double *y, double *qraux, DIMS dm, LENGTHS glen)
{ /* return the pre-decomposition for ZXy */
  QRStruct qr;
  double *rsp;
  int info = 0;

  for (int i = 0; i < dm->n; i++) {
    double *mat = ZX + (glen->offsets)[i];
    double *tau = qraux + i * dm->ZXcols;
    qr  = QR_decomp(mat, dm->ZXrows, (glen->glen)[i], dm->ZXcols, tau);
    rsp = (double *) Calloc((glen->glen)[i], double);
    Memcpy(rsp, y + (glen->offsets)[i], (glen->glen)[i]);
    QR_qty(qr, y + (glen->offsets)[i], rsp, (glen->glen)[i], 1, &info);
    if (info)
      error("QR_qty in pre_decomp gave code %d", info);
    Free(rsp); QR_free(qr);
  }
}

static int
lme_iterate(LME model)
{ /* EM estimation for lme under heavy-tailed distributions */
  int iter = 0, cycle;
  double conv, RSS, newRSS;

  /* initialization */
  RSS = (double) model->dm->N;

  /* main loop */
  repeat {
    /* internal EM cycles */
    for (cycle = 1; cycle <= model->ncycles; cycle++)
      internal_EMcycle(model);
    /* outer E-step */
    outer_Estep(model);
    newRSS  = sum_abs(model->distances, 1, model->dm->n);
    newRSS *= (model->scale)[0];

    iter++;

    /* eval convergence */
    conv = fabs((newRSS - RSS) / (newRSS + ABSTOL));
    if (conv < model->tolerance) /* successful completion */
      break;
    if (iter >= model->maxIter)
      break; /* maximum number of iterations exceeded */
    RSS = newRSS;
  }
  return iter;
}

static void
internal_EMcycle(LME model)
{
  internal_Estep(model);
  update_coef(model);
  update_theta(model);
  update_scale(model);
  if (!model->fixShape)
    update_nu(model);
}

static void
internal_Estep(LME model)
{
  DIMS dm = model->dm;
  LENGTHS glen = model->lengths;
  double *Root;

  Root = (double *) Calloc(SQR(dm->ZXcols + 1), double);
  relative_precision(model->theta, dm->q, model->scale, model->Delta);
  for (int i = 0; i < dm->n; i++) {
    double *R = model->ZX + (glen->offsets)[i];
    double *c = model->y + (glen->offsets)[i];
    append_decomp(R, c, (glen->glen)[i], dm->ZXrows, dm->ZXcols, model->Delta, dm->q, Root, dm->ZXcols + 1);
    double *ranef = model->ranef + i * dm->q;
    random_effects(Root, dm->ZXcols + 1, dm, model->coef, ranef);
    (model->distances)[i] = mahalanobis(Root, dm->ZXcols + 1, dm, model->coef, model->scale);
    double *Half = model->Root + i * dm->q * dm->q;
    upper_tri(Half, dm->q, Root, dm->ZXcols + 1, dm->q, dm->q);
    setzero(Root, dm->ZXcols + 1, dm->ZXcols + 1, dm->ZXcols + 1);
  }
  Free(Root);
}

static void
outer_Estep(LME model)
{
  DIMS dm = model->dm;
  LENGTHS glen = model->lengths;

  for (int i = 0; i < dm->n; i++)
    (model->weights)[i] = do_weight(model->family, (glen->glen)[i], (model->distances)[i]);
}

static void
update_coef(LME model)
{ /* compute the fixed effect estimates */
  DIMS dm = model->dm;
  LENGTHS glen = model->lengths;
  double wts, *R, *X, *working;
  int info = 0;

  X = (double *) Calloc(dm->DcRows * dm->p, double);
  working = (double *) Calloc(dm->DcRows, double);
  for (int i = 0; i < dm->n; i++) {
    R = (double *) Calloc((glen->DcLen)[i] * dm->ZXcols, double);
    double *ZX = model->ZX + (glen->offsets)[i];
    double *rsp = model->y + (glen->offsets)[i];
    double *ranef = model->ranef + i * dm->q;
    wts = sqrt((model->weights)[i]);
    upper_tri(R, (glen->DcLen)[i], ZX, dm->ZXrows, (glen->DcLen)[i], dm->ZXcols);
    working_response(R, rsp, (glen->DcLen)[i], dm, ranef, working + (glen->DcOff)[i]);
    scale_mat(X + (glen->DcOff)[i], dm->DcRows, wts, R + (glen->DcLen)[i] * dm->q,
              (glen->DcLen)[i], (glen->DcLen)[i], dm->p);
    scale_mat(working + (glen->DcOff)[i], 0, wts, working + (glen->DcOff)[i], 0,
              (glen->DcLen)[i], 0); // 1
    Free(R);
  }
  lsfit(X, dm->DcRows, dm->DcRows, dm->p, working, dm->DcRows, 1, model->coef, &info);
  if (info)
    error("lsfit in update_coef gave code %d", info);
  Free(X); Free(working);
}

static void
update_scale(LME model)
{ /* update the within-groups scale parameter */
  DIMS dm = model->dm;
  double SSQ;

  SSQ  = dot_product(model->distances, 1, model->weights, 1, dm->n);
  SSQ += dm->n * dm->q;
  (model->scale)[0] *= SSQ;
  (model->scale)[0] /= dm->N + dm->n * dm->q;
}

static void
update_theta(LME model)
{ /* update the scale matrix of random effects */
  DIMS dm = model->dm;
  int info = 0, job = 1, qsq = SQR(dm->q);
  double *accum, *Omega, scale, wts;

  scale = (model->scale)[0];
  accum = (double *) Calloc(qsq, double);
  Omega = (double *) Calloc(qsq, double);
  for (int i = 0; i < dm->n; i++) {
    double *Root = model->Root + i * qsq;
    invert_triangular(Root, dm->q, dm->q, job, &info);
    if (info)
      error("invert_triangular in update_theta gave code %d", info);
    outerprod(Omega, Root, dm->q, dm->q, dm->q, Root, dm->q, dm->q, dm->q);
    scale_mat(Omega, dm->q, scale, Omega, dm->q, dm->q, dm->q);
    wts = (model->weights)[i];
    double *ranef = model->ranef + i * dm->q;
    rank1_update(Omega, dm->q, dm->q, dm->q, wts, ranef, ranef);
    add_mat(accum, dm->q, 1.0, Omega, dm->q, dm->q, dm->q);
    setzero(Omega, dm->q, dm->q, dm->q);
  }
  scale_mat(model->theta, dm->q, 1.0 / dm->n, accum, dm->q, dm->q, dm->q);
  Free(accum); Free(Omega);
}

static void
update_nu(LME model)
{
  DIMS dm = model->dm;
  double *lengths;

  lengths = (double *) Calloc(dm->n, double);
  for (int i = 0; i < dm->n; i++)
    lengths[i] = (double) (model->lengths->glen)[i];
  update_mixture(model->family, dm, model->distances, lengths, model->weights, model->tolerance);
  Free(lengths);
}

static void
relative_precision(double *Psi, int q, double *scale, double *Delta)
{
  int info = 0, job = 0; /* Delta is lower triangular */

  setzero(Delta, q, q, q);
  lower_tri(Delta, q, Psi, q, q, q);
  chol_decomp(Delta, q, q, job, &info);
  if (info)
    error("chol_decomp in relative_precision gave code %d", info);
  invert_triangular(Delta, q, q, job, &info);
  if (info)
    error("invert_triangular in relative_precision gave code %d", info);
  scale_mat(Delta, q, sqrt(*scale), Delta, q, q, q);
}

static void
append_decomp(double *R, double *c, int glen, int nrow, int ncol, double *Delta,
  int q, double *Store, int ldStr)
{ /* apply a QR decomposition to the augmented matrix rbind(Rc, Delta),
   * triangular part is returned in Store */
  int arow = glen + q, acol = ncol + 1;
  double *aug, *qraux;
  QRStruct qr;

  aug   = (double *) Calloc(arow * acol, double);
  qraux = (double *) Calloc(acol, double);

  upper_tri(aug, arow, R, nrow, glen, ncol);
  Memcpy(aug + arow * ncol, c, glen);
  lower_tri(aug + glen, arow, Delta, q, q, q);

  qr = QR_decomp(aug, arow, arow, acol, qraux);
  QR_store_R(qr, Store, ldStr);

  QR_free(qr); Free(aug); Free(qraux);
}

static void
random_effects(double *Root, int ldRoot, DIMS dm, double *coef, double *ranef)
{ /* compute random effects estimates */
  int info = 0;

  Memcpy(ranef, Root + ldRoot * dm->ZXcols, dm->q);
  GE_axpy(ranef, -1.0, Root + ldRoot * dm->q, ldRoot, dm->q, dm->p, coef, 1.0, 0);
  backsolve(Root, ldRoot, dm->q, ranef, dm->q, 1, 1, &info);
  if (info)
    error("backsolve in random_effects gave code %d", info);
}

static double
mahalanobis(double *Root, int ldRoot, DIMS dm, double *coef, double *scale)
{ /* Mahalanobis distances */
  double dist, *z;

  z = (double *) Calloc(dm->p + 1, double);
  Memcpy(z, Root + ldRoot * dm->ZXcols + dm->q, dm->p + 1);
  GE_axpy(z, -1.0, Root + ldRoot * dm->q + dm->q, ldRoot, dm->p, dm->p, coef, 1.0, 0);
  dist = norm_sqr(z, 1, dm->p + 1) / *scale;
  Free(z);
  return dist;
}

static void
working_response(double *R, double *c, int DcLen, DIMS dm, double *ranef, double *working)
{ /* compute the working response */
  Memcpy(working, c, DcLen);
  GE_axpy(working, -1.0, R, DcLen, dm->q, dm->q, ranef, 1.0, 0);
}

static double
lme_logLik(LME model)
{ /* evaluate the log-likelihood function */
  DIMS dm = model->dm;
  int qsq = SQR(dm->q);
  double *lengths, accum = 0.0, kernel, logDet, scale;

  scale = (model->scale)[0];
  relative_precision(model->theta, dm->q, model->scale, model->Delta);
  lengths = (double *) Calloc(dm->n, double);
  for (int i = 0; i < dm->n; i++) {
    double *Root = model->Root + i * qsq;
    accum += logAbsDet(Root, dm->q, dm->q);
    lengths[i] = (double) (model->lengths->glen)[i];
  }
  logDet = logAbsDet(model->Delta, dm->q, dm->q);
  kernel = logLik_kernel(model->family, dm, lengths, model->distances);
  Free(lengths);
  return (kernel - 0.5 * dm->N * log(scale) + dm->n * logDet + accum);
}

void
lme_fitted(double *ZX, int *pdims, int *lengths, int *DcLen, double *coef, double *ranef,
  double *conditional, double *marginal)
{
  FITTED ans;

  ans = lme_fitted_init(ZX, pdims, lengths, DcLen, coef, ranef, conditional, marginal);
  lme_fitted_values(ans);
  lme_fitted_free(ans);
}

static FITTED
lme_fitted_init(double *ZX, int *pdims, int *lengths, int *DcLen, double *coef, double *ranef,
  double *conditional, double *marginal)
{ /* constructor for a fitted object */
  FITTED ans;

  ans = (FITTED) Calloc(1, FITTED_struct);
  ans->ZX = ZX;
  ans->dm = dims(pdims);
  ans->lengths = setLengths(ans->dm, lengths, DcLen);
  ans->coef = coef;
  ans->ranef = ranef;
  ans->conditional = conditional;
  ans->marginal = marginal;
  return(ans);
}

static void
lme_fitted_free(FITTED this)
{ /* destructor for a fitted object */
  dims_free(this->dm);
  lengths_free(this->lengths);
  Free(this);
}

static void
lme_fitted_values(FITTED object)
{
  DIMS dm = object->dm;
  LENGTHS glen = object->lengths;
  int job = 0;

  /* marginal fitted values */
  GE_axpy(object->marginal, 1.0, object->ZX + dm->ZXrows * dm->q, dm->ZXrows, dm->ZXrows, dm->p, object->coef, 0.0, job);
  /* conditional fitted values */
  Memcpy(object->conditional, object->marginal, dm->ZXrows);
  for (int i = 0; i < dm->n; i++) {
    double *Z = object->ZX + (glen->offsets)[i];
    double *yFit = object->conditional + (glen->offsets)[i];
    double *ranef = object->ranef + i * dm->q;
    GE_axpy(yFit, 1.0, Z, dm->ZXrows, (glen->glen)[i], dm->q, ranef, 1.0, job);
  }
}

void
lme_acov(double *ZX, int *pdims, int *lengths, int *DcLen, double *settings, double *Root,
  double *scale, double *control, double *acov)
{
  ACOV ans;

  ans = lme_acov_init(ZX, pdims, lengths, DcLen, settings, Root, scale, control, acov);
  lme_acov_coef(ans);
  lme_acov_free(ans);
}

static ACOV
lme_acov_init(double *ZX, int *pdims, int *lengths, int *DcLen, double *settings,
  double *Root, double *scale, double *control, double *acov)
{ /* constructor for a covariance object */
  ACOV ans;

  ans = (ACOV) Calloc(1, ACOV_struct);
  ans->ZX = ZX;
  ans->dm = dims(pdims);
  ans->lengths = setLengths(ans->dm, lengths, DcLen);
  ans->family = family_init(settings);
  ans->Root = Root;
  ans->scale = scale;
  ans->control = control;
  ans->ndraws = (int) control[3];
  ans->acov = acov;
  return(ans);
}

static void
lme_acov_free(ACOV this)
{ /* destructor for a fitted object */
  dims_free(this->dm);
  lengths_free(this->lengths);
  family_free(this->family);
  Free(this);
}

static void
lme_acov_coef(ACOV object)
{ /* evaluate the Fisher information matrix */
  DIMS dm = object->dm;
  LENGTHS glen = object->lengths;
  FAMILY family = object->family;
  int qsq = SQR(dm->q);
  double factor, *accum, *cross, *prod, *outer, *Z, *R;

  accum = (double *) Calloc(SQR(dm->p), double);
  cross = (double *) Calloc(SQR(dm->p), double);
  prod  = (double *) Calloc(dm->p * dm->q, double);
  outer = (double *) Calloc(SQR(dm->p), double);
  for (int i = 0; i < dm->n; i++) {
    R = (double *) Calloc((glen->glen)[i] * dm->ZXcols, double);
    double *ZX = object->ZX + (glen->offsets)[i];
    upper_tri(R, (glen->glen)[i], ZX, dm->ZXrows, (glen->glen)[i], dm->ZXcols);
    crossprod(cross, R + (glen->glen)[i] * dm->q, (glen->glen)[i], (glen->glen)[i],
              dm->p, R + (glen->glen)[i] * dm->q, (glen->glen)[i], (glen->glen)[i],
              dm->p);
    Z = (double *) Calloc((glen->glen)[i] * dm->q, double);
    double *Root = object->Root + i * qsq;
    mult_mat(Z, R, (glen->glen)[i], (glen->glen)[i], dm->q, Root, dm->q, dm->q, dm->q);
    crossprod(prod, R + (glen->glen)[i] * dm->q, (glen->glen)[i], dm->q, dm->p,
              Z, (glen->glen)[i], dm->q, dm->q);
    outerprod(outer, prod, dm->p, dm->p, dm->q, prod, dm->p, dm->p, dm->q);
    add_mat(cross, dm->p, -1.0, outer, dm->p, dm->p, dm->p);
    factor  = acov_scale(family, (glen->glen)[i], object->ndraws);
    factor /= (glen->glen)[i];
    add_mat(accum, dm->p, factor, cross, dm->p, dm->p, dm->p);
    Free(Z); Free(R);
  }
  copy_mat(object->acov, dm->p, accum, dm->p, dm->p, dm->p);
  Free(accum); Free(prod); Free(cross); Free(outer);
}
