#include <R.h>
#include <Rinternals.h>

#include <R_ext/Rdynload.h>

SEXP fgsea_calcGseaStatCumulative(SEXP statsSEXP, SEXP nSEXP, SEXP kSEXP, SEXP gseaParamSEXP);

R_CallMethodDef callMethods[]  = {
  {"fgsea_calcGseaStatCumulative", (DL_FUNC) &fgsea_calcGseaStatCumulative, 4},
  {NULL, NULL, 0}
};

void R_init_fgsea(DllInfo *info) {
  R_registerRoutines(info, NULL, callMethods, NULL, NULL);
  R_useDynamicSymbols(info, FALSE);
}

