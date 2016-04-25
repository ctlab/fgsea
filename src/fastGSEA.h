#include <Rcpp.h>
using namespace Rcpp;

#include <vector>
using namespace std;


// [[Rcpp::export]]
NumericVector calcGseaStatCumulative(
        NumericVector const& stats, 
        IntegerVector const& selectedStats, // Indexes start from one!
        double gseaParam
        );
