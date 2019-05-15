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

    EsRuler esRulerPos(posRanks, sampleSize, pathwaySize);
    EsRuler esRulerNeg(negRanks, sampleSize, pathwaySize);

    double maxES = *max_element(begin(esVector), end(esVector));
    double minES = *min_element(begin(esVector), end(esVector));
    if (maxES >= 0){
        esRulerPos.extend(abs(maxES), seed, absEps);
    }
    if (minES < 0){
        esRulerNeg.extend(abs(minES), seed, absEps);
    }

    vector<double> result;
    int nrow = esVector.size();
    for (int i = 0; i < nrow; i++){
        double pvalue;
        double currentES = esVector[i];
        if (currentES >= 0.0){
            pvalue = esRulerPos.getPvalue(abs(currentES), absEps, sign);
            result.emplace_back(pvalue);
        }
        else{
            pvalue = esRulerNeg.getPvalue(abs(currentES), absEps, sign);
            result.emplace_back(pvalue);
        }
    }
    return wrap(result);
}
