#include "fgseaMultilevel.h"
#include "fgseaMultilevelSupplement.h"

#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
using namespace std;


NumericVector fgseaMultilevelCpp(const NumericVector& inpSizes, const NumericVector inpEs, 
                               const NumericVector& ranks, int samplesSize,
                               int seed, double absEps, bool sign)
{
    vector<double> result;

    vector<double> S_pos = as<std::vector<double> >(ranks);
    for (int i = 0; i < S_pos.size(); i++) {
        S_pos[i] = fabs(S_pos[i]);
    }
    vector<double> S_neg = S_pos;
    reverse(S_neg.begin(), S_neg.end());

    const vector<double> esVector = as<std::vector<double> >(inpEs);
    const vector<int> sizesVector = as<std::vector<int> >(inpSizes);

    int nrow = esVector.size();
    int i = 0;
    while(i < nrow)
    {
        double pvalue;
        double ES = esVector[i];
        int pathwaySize = sizesVector[i];
        vector<double> &S = (ES > 0 ? S_pos : S_neg);
        EsPvalConnection esPvalObj(samplesSize);

        calcPvalues(esPvalObj, S, pathwaySize, fabs(ES), samplesSize, seed, absEps);
        if (ES < 0)
        {
            while (esVector[i] < 0 && (i < nrow) && (sizesVector[i] == pathwaySize))
            {
                pvalue = findEsPval(esPvalObj, fabs(esVector[i]), samplesSize, sign);
                result.emplace_back(pvalue);
                i++;
            }
        }
        else
        {
            while (esVector[i] > 0 && (i < nrow) && (sizesVector[i] == pathwaySize))
            {
                pvalue = findEsPval(esPvalObj, fabs(esVector[i]), samplesSize, sign);
                result.emplace_back(pvalue);
                i++;
            }
        }
    }
    return wrap(result);
}
