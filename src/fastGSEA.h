#include <Rcpp.h>
using namespace Rcpp;

#include <vector>
using namespace std;


//' Calculates GSEA statistic values for all the prefixes of a gene set
//' @param stats Named numeric vector with gene-level statistics
//' @param selectedStats indexes of selected genes in a 'stats' array
//' @param gseaParam GSEA weight parameter (0 is unweighted, suggested value is 1)
//' @export
// [[Rcpp::export]]
NumericVector calcGseaStatCumulative(
        NumericVector const& stats,
        IntegerVector const& selectedStats, // Indexes start from one!
        double gseaParam
        );
