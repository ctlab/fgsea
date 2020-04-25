#include <R.h>
#include <Rinternals.h>

#include <R_ext/Rdynload.h>

SEXP _fgsea_calcGseaStatCumulativeBatch(SEXP statsSEXP, SEXP gseaParamSEXP, SEXP pathwayScoresSEXP, SEXP pathwaysSizesSEXP, SEXP iterationsSEXP, SEXP seedSEXP, SEXP scoreTypeSEXP);
SEXP _fgsea_calcGseaStatCumulative(SEXP statsSEXP, SEXP selectedStatsSEXP, SEXP gseaParamSEXP, SEXP scoreTypeSEXP);
SEXP _fgsea_calcGseaStatBatchCpp(SEXP statsSEXP, SEXP selectedGenesSEXP, SEXP geneRanksSEXP);
SEXP _fgsea_fgseaMultilevelCpp(SEXP enrichmentScoresSEXP, SEXP ranksSEXP, SEXP pathwaySizeSEXP, SEXP sampleSizeSEXP, SEXP seedSEXP, SEXP epsSEXP, SEXP signSEXP);

R_CallMethodDef callMethods[]  = {
  {"_fgsea_calcGseaStatCumulativeBatch", (DL_FUNC) &_fgsea_calcGseaStatCumulativeBatch, 7},
  {"_fgsea_calcGseaStatCumulative", (DL_FUNC) &_fgsea_calcGseaStatCumulative, 4},
  {"_fgsea_calcGseaStatBatchCpp", (DL_FUNC) &_fgsea_calcGseaStatBatchCpp, 3},
  {"_fgsea_fgseaMultilevelCpp", (DL_FUNC) &_fgsea_fgseaMultilevelCpp, 7},
  {NULL, NULL, 0}
};

void R_init_fgsea(DllInfo *info) {
  R_registerRoutines(info, NULL, callMethods, NULL, NULL);
  R_useDynamicSymbols(info, FALSE);
}

