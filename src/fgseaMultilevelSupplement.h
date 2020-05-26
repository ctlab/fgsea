#ifndef FGSEAMULTILEVELCPP_FGSEAMULTILEVELSUPPLEMENT_H
#define FGSEAMULTILEVELCPP_FGSEAMULTILEVELSUPPLEMENT_H

#include <vector>
#include <random>
#include <set>
#include <cmath>
#include <algorithm>
#include <boost/math/special_functions/digamma.hpp>

#include <iostream>

#include <ext/pb_ds/assoc_container.hpp>

using namespace __gnu_pbds;
using namespace std;

template <typename T> using ordered_set = tree<T, null_type, less<T>, rb_tree_tag, tree_order_statistics_node_update>;
template <typename K, typename V> using ordered_map = tree<K, V, less<K>, rb_tree_tag, tree_order_statistics_node_update>;


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

    vector<int> chunkLastElement;
    int chunksNumber;

    struct SampleChunks {
        vector<double> chunkSum;
        vector<int> chunkSize;
        vector<vector<int>> chunks;
        SampleChunks(int);
    };

    int perturbate(const vector<double> &ranks, int k, SampleChunks &cusSampleChunks,
               double bound, mt19937 &rng);

    int chunkLen(int ind);

public:

    EsRuler(const vector<double> &inpRanks, unsigned int inpSampleSize, unsigned int inpPathwaySize);

    ~EsRuler();

    void extend(double ES, int seed, double eps);

    pair<double, bool> getPvalue(double ES, double eps, bool sign);
};

double betaMeanLog(unsigned long a, unsigned long b);

pair<double, bool> calcLogCorrection(const vector<unsigned int> &probCorrector, long probCorrIndx,
                         const pair<unsigned int, unsigned int> posUnifScoreCount, unsigned int sampleSize);


#endif //FGSEAMULTILEVELCPP_FGSEAMULTILEVELSUPPLEMENT_H
