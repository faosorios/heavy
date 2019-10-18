/* ID: matrix.h, last updated 2019/06/07, F. Osorio */

#ifndef HEAVY_MATRIX_H
#define HEAVY_MATRIX_H

#include "base.h"

/* basic vector manipulations and BLAS-1 wrappers */
extern double dot_product(double *, int, double *, int, int);
extern double norm_sqr(double *, int, int);
extern double sum_abs(double *, int, int);
extern void ax_plus_y(double, double *, int, double *, int, int);
extern void copy_vec(double *, int, double *, int, int);
extern void scale_vec(double, double *, int, int);
extern void swap_vec(double *, int, double *, int, int);
extern void hadamard_prod(double *, double *, double *, int);

/* basic matrix manipulations and BLAS-2, and BLAS-3 wrappers */
extern void add_mat(double *, int, double, double *, int, int, int);
extern void copy_mat(double *, int, double *, int, int, int);
extern void crossprod(double *, double *, int, int, int, double *, int, int, int);
extern void GE_axpy(double *, double, double *, int, int, int, double *, double, int);
extern void lower_tri(double *, int, double *, int, int, int);
extern void mult_mat(double *, double *, int, int, int, double *, int, int, int);
extern void outerprod(double *, double *, int, int, int, double *, int, int, int);
extern void rank1_update(double *, int, int, int, double, double *, double *);
extern void rank1_symm_update(double *, int, int, double, double *, int);
extern void scale_mat(double *, int, double, double *, int, int, int);
extern void setzero(double *, int, int, int);
extern void triangle_mult_vec(double *, double *, int, int, double *, int);
extern void upper_tri(double *, int, double *, int, int, int);
extern double logAbsDet(double *, int, int);

/* DEBUG routine */
extern void print_mat(double *, int, int, int, char *);

/* routines for matrix decompositions (wrappers to LAPACK and Linpack) */
extern void chol_decomp(double *, int, int, int, int *);
extern void svd_decomp(double *, int, int, int, double *, int, double *, double *, int, int, int *);
extern QRStruct QR_decomp(double *, int, int, int, double *);
extern void QR_free(QRStruct);

/* orthogonal-triangular operations (wrappers to Linpack) */
extern void QR_qty(QRStruct, double *, double *, int, int, int *);
extern void QR_qy(QRStruct, double *, double *, int, int, int *);
extern void QR_coef(QRStruct, double *, double *, int, int, int *);
extern void QR_resid(QRStruct, double *, double *, int, int, int *);
extern void QR_fitted(QRStruct, double *, double *, int, int, int *);
extern void QR_store_R(QRStruct, double *, int);

/* matrix inversion and linear solver */
extern void invert_mat(double *, int, int, int *);
extern void invert_triangular(double *, int, int, int, int *);
extern void backsolve(double *, int, int, double *, int, int, int, int *);

/* linear least-squares fit */
extern void lsfit(double *, int, int, int, double *, int, int, double *, int *);

#endif /* HEAVY_MATRIX_H */
