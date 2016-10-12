#include <Rcpp.h>
using namespace Rcpp;

#include <vector>
using namespace std;

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
List calcGseaStatCumulativeBatch(
        NumericVector const& stats,
        int n,
        int k,
        double gseaParam,
        int m,
        NumericVector const& pathwayScores,
        IntegerVector const& pathwaysSizes,
        int iterations,
        int seed);
