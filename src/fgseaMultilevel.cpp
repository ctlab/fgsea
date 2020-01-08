#include "fgseaMultilevel.h"
#include "fgseaMultilevelSupplement.h"
using namespace std;

DataFrame fgseaMultilevelCpp(const NumericVector& enrichmentScores,
                             const NumericVector& ranks, int pathwaySize,
                             int sampleSize, int seed,double eps, bool sign)
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
        esRulerPos.extend(abs(maxES), seed, eps);
    }
    if (minES < 0){
        esRulerNeg.extend(abs(minES), seed, eps);
    }

    vector<double> pvalRes;
    vector<bool> isCpGeHalf;

    int nrow = esVector.size();
    for (int i = 0; i < nrow; i++){
        pair<double, bool> resPair;
        double currentES = esVector[i];
        if (currentES >= 0.0){
            resPair = esRulerPos.getPvalue(abs(currentES), eps, sign);
            pvalRes.push_back(resPair.first);
            isCpGeHalf.push_back(resPair.second);
        }
        else{
            resPair = esRulerNeg.getPvalue(abs(currentES), eps, sign);
            pvalRes.push_back(resPair.first);
            isCpGeHalf.push_back(resPair.second);
        }
    }

    // return vector with pvalues and vector with conditional probability result
    return DataFrame::create(Named("cppMPval") = pvalRes,
                             Named("cppIsCpGeHalf") = isCpGeHalf);
}
