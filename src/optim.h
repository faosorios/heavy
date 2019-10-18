/* ID: optim.h, last updated 2019/08/02, F. Osorio */

#ifndef HEAVY_OPTIM_H
#define HEAVY_OPTIM_H

#include "base.h"

/* Brent's method for unidimensional optimization */
extern double brent(double, double, double (*f)(double, void *), void *, double);

#endif /* HEAVY_OPTIM_H */
