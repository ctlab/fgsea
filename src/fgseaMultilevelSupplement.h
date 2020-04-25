#ifndef FGSEAMULTILEVELCPP_FGSEAMULTILEVELSUPPLEMENT_H
#define FGSEAMULTILEVELCPP_FGSEAMULTILEVELSUPPLEMENT_H

#include <vector>
#include <random>
#include <set>
#include <cmath>
#include <algorithm>
#include <boost/math/special_functions/digamma.hpp>

using namespace std;


class EsRuler {
private:

    const vector<double> &ranks;
    const unsigned int sampleSize;
    const unsigned int pathwaySize;

    //first for positive ES count
    //second for total count
    pair<unsigned int, unsigned int> posUnifScoreCount;
    
    vector<double> enrichmentScores;
    vector<vector<int> > currentSamples;

    vector<unsigned int> probCorrector;

    void duplicateSamples();

public:

    EsRuler(const vector<double> &inpRanks, unsigned int inpSampleSize, unsigned int inpPathwaySize);

    ~EsRuler();

    void extend(double ES, int seed, double eps);

    pair<double, bool> getPvalue(double ES, double eps, bool sign);
};

int perturbate(const vector<double> &ranks, vector<int> &sample,
               double bound, mt19937 &rng);

double betaMeanLog(unsigned long a, unsigned long b);

pair<double, bool> calcLogCorrection(const vector<unsigned int> &probCorrector, long probCorrIndx,
                         const pair<unsigned int, unsigned int> posUnifScoreCount, unsigned int sampleSize);


#endif //FGSEAMULTILEVELCPP_FGSEAMULTILEVELSUPPLEMENT_H
