#include <Rcpp.h>
using namespace Rcpp;

#include <vector>
using namespace std;

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
List calcGseaStatCumulativeBatch(
        NumericVector const& stats,
        double gseaParam,
        NumericVector const& pathwayScores,
        IntegerVector const& pathwaysSizes,
        int iterations,
        int seed);

//' Calculates GSEA statistic values for all the prefixes of a gene set
//'
//' Takes \emph{O(k^\{3/2\})} time, where \emph{k} is a size of `selectedSize`.
//' @param stats Named numeric vector with gene-level statistics
//'     sorted in decreasing order (order is not checked)
//' @param selectedStats indexes of selected genes in a `stats` array
//' @param gseaParam GSEA weight parameter (0 is unweighted, suggested value is 1)
//' @return Numeric vector of GSEA statistics for all prefixes of selectedStats.
//' @export
//' @examples
//' data(exampleRanks)
//' data(examplePathways)
//' ranks <- sort(exampleRanks, decreasing=TRUE)
//' es <- calcGseaStatCumulative(ranks, na.omit(match(examplePathways[[1]], names(ranks))), 1)
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
NumericVector calcGseaStatCumulative(
        NumericVector const& stats,
        IntegerVector const& selectedStats, // Indexes start from one!
        double gseaParam
);
