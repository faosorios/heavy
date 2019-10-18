/* ID: random.c, last updated 2019/08/02, F. Osorio */

#include "random.h"

/* static functions.. */
static DIMS dims(int *);
static void dims_free(DIMS);
/* ..end declarations */

/* 'dims' functions */

static DIMS
dims(int *pdims)
{ /* dims object */
  DIMS ans;

  ans = (DIMS) Calloc(1, DIMS_struct);
  ans->n = (int) pdims[0];
  ans->p = (int) pdims[1];
  return ans;
}

static void
dims_free(DIMS this)
{ /* destructor for a dims object */
  Free(this);
}

/* random vector generation uniformly located on a spherical surface */

void
rand_sphere(double *y, int *pdims)
{ /* random vector generation on the sphere (to be called by R) */
  DIMS dm;

  dm = dims(pdims);
  GetRNGstate();
  rand_unif_sphere(y, dm->n, dm->p);
  PutRNGstate();
  dims_free(dm);
}

void
rand_unif_sphere(double *y, int n, int p)
{ /* random vector generation uniformly on the sphere */
  int inc = 1;
  double radial;

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < p; j++)
      y[j] = norm_rand();
    radial = 1.0 / F77_CALL(dnrm2)(&p, y, &inc);
    F77_CALL(dscal)(&p, &radial, y, &inc);
    y += p;
  }
}

/* multivariate normal random generation */

void
rand_norm(double *y, int *pdims, double *center, double *Scatter)
{ /* multivariate normal random generation */
  DIMS dm;
  char *side = "L", *uplo = "U", *trans = "T", *diag = "N";
  double one = 1.;
  int inc = 1, info = 0, job = 1;

  dm = dims(pdims);
  GetRNGstate();
  chol_decomp(Scatter, dm->p, dm->p, job, &info);
  if (info)
    error("cholesky factorization in rand_norm gave code %d", info);
  rand_spherical_norm(y, dm->n, dm->p);
  F77_CALL(dtrmm)(side, uplo, trans, diag, &(dm->p), &(dm->n), &one, Scatter, &(dm->p), y, &(dm->p));
  for (int i = 0; i < dm->n; i++) {
    F77_CALL(daxpy)(&(dm->p), &one, center, &inc, y, &inc);
    y += dm->p;
  }
  PutRNGstate();
  dims_free(dm);
}

void
rand_spherical_norm(double *y, int n, int p)
{ /* independent standard normal variates */
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < p; j++)
      y[j] = norm_rand();
    y += p;
  }
}

/* multivariate Cauchy random generation */

void
rand_cauchy(double *y, int *pdims, double *center, double *Scatter)
{ /* multivariate normal random generation */
  DIMS dm;
  char *side = "L", *uplo = "U", *trans = "T", *diag = "N";
  double one = 1.;
  int inc = 1, info = 0, job = 1;

  dm = dims(pdims);
  GetRNGstate();
  chol_decomp(Scatter, dm->p, dm->p, job, &info);
  if (info)
    error("cholesky factorization in rand_cauchy gave code %d", info);
  rand_spherical_cauchy(y, dm->n, dm->p);
  F77_CALL(dtrmm)(side, uplo, trans, diag, &(dm->p), &(dm->n), &one, Scatter, &(dm->p), y, &(dm->p));
  for (int i = 0; i < dm->n; i++) {
    F77_CALL(daxpy)(&(dm->p), &one, center, &inc, y, &inc);
    y += dm->p;
  }
  PutRNGstate();
  dims_free(dm);
}

void
rand_spherical_cauchy(double *y, int n, int p)
{ /* standard Cauchy variates */
  int inc = 1;
  double tau, radial;

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < p; j++)
      y[j] = norm_rand();
    tau = rgamma(0.5, 2.0);
    radial = R_pow(tau, -0.5);
    F77_CALL(dscal)(&p, &radial, y, &inc);
    y += p;
  }
}

/* multivariate Student-t random generation */

void
rand_student(double *y, int *pdims, double *center, double *Scatter, double *df)
{ /* multivariate Student-t random generation */
  DIMS dm;
  char *side = "L", *uplo = "U", *trans = "T", *diag = "N";
  double one = 1.;
  int inc = 1, info = 0, job = 1;

  dm = dims(pdims);
  GetRNGstate();
  chol_decomp(Scatter, dm->p, dm->p, job, &info);
  if (info)
    error("cholesky factorization in rand_student gave code %d", info);
  rand_spherical_student(y, *df, dm->n, dm->p);
  F77_CALL(dtrmm)(side, uplo, trans, diag, &(dm->p), &(dm->n), &one, Scatter, &(dm->p), y, &(dm->p));
  for (int i = 0; i < dm->n; i++) {
    F77_CALL(daxpy)(&(dm->p), &one, center, &inc, y, &inc);
    y += dm->p;
  }
  PutRNGstate();
  dims_free(dm);
}

void
rand_spherical_student(double *y, double df, int n, int p)
{ /* standard Student-t variates */
  int inc = 1;
  double tau, radial;

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < p; j++)
      y[j] = norm_rand();
    tau = rgamma(df / 2., 2. / df);
    radial = R_pow(tau, -.5);
    F77_CALL(dscal)(&p, &radial, y, &inc);
    y += p;
  }
}

/* multivariate slash random generation */

void
rand_slash(double *y, int *pdims, double *center, double *Scatter, double *df)
{ /* multivariate slash random generation */
  DIMS dm;
  char *side = "L", *uplo = "U", *trans = "T", *diag = "N";
  double one = 1.;
  int inc = 1, info = 0, job = 1;

  dm = dims(pdims);
  GetRNGstate();
  chol_decomp(Scatter, dm->p, dm->p, job, &info);
  if (info)
    error("cholesky factorization in rand_slash gave code %d", info);
  rand_spherical_slash(y, *df, dm->n, dm->p);
  F77_CALL(dtrmm)(side, uplo, trans, diag, &(dm->p), &(dm->n), &one, Scatter, &(dm->p), y, &(dm->p));
  for (int i = 0; i < dm->n; i++) {
    F77_CALL(daxpy)(&(dm->p), &one, center, &inc, y, &inc);
    y += dm->p;
  }
  PutRNGstate();
  dims_free(dm);
}

void
rand_spherical_slash(double *y, double df, int n, int p)
{ /* standard slash variates */
  int inc = 1;
  double tau, radial;

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < p; j++)
      y[j] = norm_rand();
    tau = rbeta(df, 1.);
    radial = R_pow(tau, -.5);
    F77_CALL(dscal)(&p, &radial, y, &inc);
    y += p;
  }
}

/* multivariate contaminated normal random generation */

void
rand_contaminated(double *y, int *pdims, double *center, double *Scatter, double *eps, double *vif)
{ /* multivariate slash random generation */
  DIMS dm;
  char *side = "L", *uplo = "U", *trans = "T", *diag = "N";
  double one = 1.;
  int inc = 1, info = 0, job = 1;

  dm = dims(pdims);
  GetRNGstate();
  chol_decomp(Scatter, dm->p, dm->p, job, &info);
  if (info)
    error("cholesky factorization in rand_contaminated gave code %d", info);
  rand_spherical_contaminated(y, *eps, *vif, dm->n, dm->p);
  F77_CALL(dtrmm)(side, uplo, trans, diag, &(dm->p), &(dm->n), &one, Scatter, &(dm->p), y, &(dm->p));
  for (int i = 0; i < dm->n; i++) {
    F77_CALL(daxpy)(&(dm->p), &one, center, &inc, y, &inc);
    y += dm->p;
  }
  PutRNGstate();
  dims_free(dm);
}

void
rand_spherical_contaminated(double *y, double eps, double vif, int n, int p)
{ /* standard contaminated normal variates */
  int inc = 1;
  double radial, unif;

  radial = 1. / vif;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < p; j++)
      y[j] = norm_rand();
    unif = unif_rand();
    if (unif > 1. - eps)
      F77_CALL(dscal)(&p, &radial, y, &inc);
    y += p;
  }
}
