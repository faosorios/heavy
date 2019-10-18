/* ID: random.h, last updated 2019/08/02, F. Osorio */

#ifndef HEAVY_RANDOM_H
#define HEAVY_RANDOM_H

#include "matrix.h"

/* multivariate symmetric random generation (to be called by R) */
extern void rand_sphere(double *, int *);
extern void rand_norm(double *, int *, double *, double *);
extern void rand_cauchy(double *, int *, double *, double *);
extern void rand_student(double *, int *, double *, double *, double *);
extern void rand_slash(double *, int *, double *, double *, double *);
extern void rand_contaminated(double *, int *, double *, double *, double *, double *);

/* spherical random generation */
extern void rand_spherical_norm(double *, int, int);
extern void rand_spherical_cauchy(double *, int, int);
extern void rand_spherical_student(double *, double, int, int);
extern void rand_spherical_slash(double *, double, int, int);
extern void rand_spherical_contaminated(double *, double, double, int, int);

/* uniformly distributed random vectors on the unitary sphere */
extern void rand_unif_sphere(double *, int, int);

#endif /* HEAVY_RANDOM_H */
