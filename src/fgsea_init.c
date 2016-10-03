#include <R.h>
#include <Rinternals.h>

#include <R_ext/Rdynload.h>

SEXP fgsea_calcGseaStatCumulativeParallel(SEXP statsSEXP, SEXP nSEXP, SEXP kSEXP, SEXP gseaParamSEXP, SEXP mSEXP, SEXP pathwayScoresSEXP, SEXP pathwaysSizesSEXP, SEXP iterationsSEXP, SEXP setSeedSEXP, SEXP seedSEXP);

R_CallMethodDef callMethods[]  = {
  {"fgsea_calcGseaStatCumulativeParallel", (DL_FUNC) &fgsea_calcGseaStatCumulativeParallel, 10},
  {NULL, NULL, 0}
};

void R_init_fgsea(DllInfo *info) {
  R_registerRoutines(info, NULL, callMethods, NULL, NULL);
  R_useDynamicSymbols(info, FALSE);
}

