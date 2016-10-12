#include <R.h>
#include <Rinternals.h>

#include <R_ext/Rdynload.h>

SEXP fgsea_calcGseaStatCumulativeBatch(SEXP statsSEXP, SEXP nSEXP, SEXP kSEXP, SEXP gseaParamSEXP, SEXP mSEXP, SEXP pathwayScoresSEXP, SEXP pathwaysSizesSEXP, SEXP iterationsSEXP, SEXP seedSEXP);

R_CallMethodDef callMethods[]  = {
  {"fgsea_calcGseaStatCumulativeBatch", (DL_FUNC) &fgsea_calcGseaStatCumulativeBatch, 9},
  {NULL, NULL, 0}
};

void R_init_fgsea(DllInfo *info) {
  R_registerRoutines(info, NULL, callMethods, NULL, NULL);
  R_useDynamicSymbols(info, FALSE);
}

