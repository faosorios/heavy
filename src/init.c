/* ID: init.c, last updated 2019/08/02, F. Osorio */

#include <R_ext/Rdynload.h>
#include "distn.h"
#include "grubbs_fit.h"
#include "lm_fit.h"
#include "lme_fit.h"
#include "mv_fit.h"
#include "ps_fit.h"
#include "matrix.h" /* DEBUG */
#include "random.h"
#include "specfun.h"

static const R_CMethodDef CEntries[]  = {
  {"pdf_tgamma",            (DL_FUNC) &pdf_tgamma,            10},
  {"cdf_tgamma",            (DL_FUNC) &cdf_tgamma,            10},
  {"quantile_tgamma",       (DL_FUNC) &quantile_tgamma,       10},
  {"rand_tgamma",           (DL_FUNC) &rand_tgamma,            8},
  {"cdf_gamma_derivatives", (DL_FUNC) &cdf_gamma_derivatives,  4},
  {"grubbs_fit",            (DL_FUNC) &grubbs_fit,            13},
  {"lm_fit",                (DL_FUNC) &lm_fit,                13},
  {"lme_fit",               (DL_FUNC) &lme_fit,               16},
  {"lme_fitted",            (DL_FUNC) &lme_fitted,             8},
  {"lme_acov",              (DL_FUNC) &lme_acov,               9},
  {"mlm_fit",               (DL_FUNC) &mlm_fit,               13},
  {"mv_fit",                (DL_FUNC) &mv_fit,                10},
  {"ps_fit",                (DL_FUNC) &ps_fit,                17},
  {"ps_combined",           (DL_FUNC) &ps_combined,           17},
  {"rand_cauchy",           (DL_FUNC) &rand_cauchy,            4},
  {"rand_contaminated",     (DL_FUNC) &rand_contaminated,      6},
  {"rand_norm",             (DL_FUNC) &rand_norm,              4},
  {"rand_slash",            (DL_FUNC) &rand_slash,             5},
  {"rand_sphere",           (DL_FUNC) &rand_sphere,            2},
  {"rand_student",          (DL_FUNC) &rand_student,           5},
  {NULL, NULL, 0}
};

void R_init_heavy(DllInfo *dll)
{
  R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
