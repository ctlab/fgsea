#pragma once 

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::export]]
NumericVector fgseaMultilevelCpp(const NumericVector& inpSizes, const NumericVector inpEs, 
                               const NumericVector& ranks, int samplesSize,
                               int seed, double absEps, bool sign);

