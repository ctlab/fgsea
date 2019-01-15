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
    /*  
    The vector `enrichmentScores` must be sorted in descending 
    order of the absolute value of the `ES`.
    */
    vector<double> posRanks = as<std::vector<double> >(ranks);
    for (int i = 0; i < posRanks.size(); i++) {
        posRanks[i] = abs(posRanks[i]);
    }
    vector<double> negRanks = posRanks;
    reverse(negRanks.begin(), negRanks.end());

    const vector<double> esVector = as<std::vector<double> >(enrichmentScores);

    double posES = 0.0;
    double negES = 0.0;
    EsPvalConnection epcPosSide(sampleSize);
    EsPvalConnection epcNegSide(sampleSize);

    vector<double> result;
    int nrow = esVector.size();
    for (int i = 0; i < nrow; i++){
        double pvalue;
        double currentES = esVector[i];
        if (currentES > posES){
            posES = currentES;
            calcPvalues(epcPosSide, posRanks, pathwaySize, abs(currentES), sampleSize, seed, absEps);
            pvalue = findEsPval(epcPosSide, abs(currentES), sampleSize, sign);
            result.emplace_back(pvalue);
        }
        else if (currentES < negES){
            negES = currentES;
            calcPvalues(epcNegSide, negRanks, pathwaySize, abs(currentES), sampleSize, seed, absEps);
            pvalue = findEsPval(epcNegSide, abs(currentES), sampleSize, sign);
            result.emplace_back(pvalue);
        }
        else if (currentES >= 0.0){
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