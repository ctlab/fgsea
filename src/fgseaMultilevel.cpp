#include "fgseaMultilevel.h"
#include "fgseaMultilevelSupplement.h"
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
using namespace std;

NumericVector fgseaMultilevelCpp(const NumericVector& enrichmentScores,
                                 const NumericVector& ranks, int pathwaySize,
                                 int sampleSize, int seed,
                                 double absEps, bool sign)
{
    vector<double> posRanks = as<std::vector<double> >(ranks);
    for (int i = 0; i < posRanks.size(); i++) {
        posRanks[i] = abs(posRanks[i]);
    }
    vector<double> negRanks = posRanks;
    reverse(negRanks.begin(), negRanks.end());

    const vector<double> esVector = as<std::vector<double> >(enrichmentScores);

    EsPvalConnection epcPosSide(sampleSize);
    EsPvalConnection epcNegSide(sampleSize);

    double maxES = *max_element(begin(esVector), end(esVector));
    double minES = *min_element(begin(esVector), end(esVector));
    if (maxES >= 0){
        calcPvalues(epcPosSide, posRanks, pathwaySize, abs(maxES), sampleSize, seed, absEps);
    }
    if (minES < 0){
        calcPvalues(epcNegSide, negRanks, pathwaySize, abs(minES), sampleSize, seed, absEps);
    }

    vector<double> result;
    int nrow = esVector.size();
    for (int i = 0; i < nrow; i++){
        double pvalue;
        double currentES = esVector[i];
        if (currentES >= 0.0){
            pvalue = findEsPval(epcPosSide, abs(currentES), sampleSize, sign);
            result.emplace_back(pvalue);
        }
        else{
            pvalue = findEsPval(epcNegSide, abs(currentES), sampleSize, sign);
            result.emplace_back(pvalue);
        }
    }
    return wrap(result);
}
