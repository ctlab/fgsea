#include <R.h>
#include <Rinternals.h>

#include <R_ext/Rdynload.h>

SEXP fgsea_calcGseaStatCumulativeBatch(SEXP statsSEXP, SEXP gseaParamSEXP, SEXP pathwayScoresSEXP, SEXP pathwaysSizesSEXP, SEXP iterationsSEXP, SEXP seedSEXP);
SEXP fgsea_calcGseaStatCumulative(SEXP statsSEXP, SEXP selectedStatsSEXP, SEXP gseaParamSEXP);

R_CallMethodDef callMethods[]  = {
  {"fgsea_calcGseaStatCumulativeBatch", (DL_FUNC) &fgsea_calcGseaStatCumulativeBatch, 6},
  {"fgsea_calcGseaStatCumulative", (DL_FUNC) &fgsea_calcGseaStatCumulative, 3},
  {NULL, NULL, 0}
};

void R_init_fgsea(DllInfo *info) {
  R_registerRoutines(info, NULL, callMethods, NULL, NULL);
  R_useDynamicSymbols(info, FALSE);
}

