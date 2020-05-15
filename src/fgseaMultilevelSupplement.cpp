#include "fgseaMultilevelSupplement.h"
#include "esCalculation.h"
#include "util.h"

#include <iostream>

double betaMeanLog(unsigned long a, unsigned long b) {
    return boost::math::digamma(a) - boost::math::digamma(b + 1);
}

pair<double, bool> calcLogCorrection(const vector<unsigned int> &probCorrector, long probCorrIndx,
                         const pair<unsigned int, unsigned int> posUnifScoreCount, unsigned int sampleSize){
    double result = 0.0;
    result -= betaMeanLog(posUnifScoreCount.first, posUnifScoreCount.second);

    unsigned long halfSize = (sampleSize + 1) / 2;
    unsigned long remainder = sampleSize - probCorrIndx % (halfSize);

    double condProb = betaMeanLog(probCorrector[probCorrIndx] + 1, remainder);
    result += condProb;

    if (exp(condProb) >= 0.5){
        return make_pair(result, true);
    }
    else{
        return make_pair(result, false);
    }
}

EsRuler::EsRuler(const vector<double> &inpRanks, unsigned int inpSampleSize, unsigned int inpPathwaySize) :
    ranks(inpRanks), sampleSize(inpSampleSize), pathwaySize(inpPathwaySize) {
    posUnifScoreCount = make_pair(0, 0);
    currentSamples.resize(inpSampleSize);
}

EsRuler::~EsRuler() = default;


void EsRuler::duplicateSamples() {
    /*
     * Removes samples with an enrichment score less than the median value and
     * replaces them with samples with an enrichment score greater than the median
     * value
    */
    vector<pair<double, int> > stats(sampleSize);
    vector<int> posEsIndxs;
    int totalPosEsCount = 0;

    for (int sampleId = 0; sampleId < sampleSize; sampleId++) {
        double sampleEsPos = calcPositiveES(ranks, currentSamples[sampleId]);
        double sampleEs = calcES(ranks, currentSamples[sampleId]);
        if (sampleEs > 0) {
            totalPosEsCount++;
            posEsIndxs.push_back(sampleId);
        }
        stats[sampleId] = make_pair(sampleEsPos, sampleId);
    }
    sort(stats.begin(), stats.end());
    for (int sampleId = 0; 2 * sampleId < sampleSize; sampleId++) {
        enrichmentScores.push_back(stats[sampleId].first);
        if (find(posEsIndxs.begin(), posEsIndxs.end(), stats[sampleId].second) != posEsIndxs.end()) {
            totalPosEsCount--;
        }
        probCorrector.push_back(totalPosEsCount);
    }

    vector<vector<int> > new_sets;
    for (int sampleId = 0; 2 * sampleId < sampleSize - 2; sampleId++) {
        for (int rep = 0; rep < 2; rep++) {
            new_sets.push_back(currentSamples[stats[sampleSize - 1 - sampleId].second]);
        }
    }
    new_sets.push_back(currentSamples[stats[sampleSize >> 1].second]);
    swap(currentSamples, new_sets);
}

void EsRuler::extend(double ES, int seed, double eps) {
    unsigned int posCount = 0;
    unsigned int totalCount = 0;
    mt19937 gen(seed);

    for (int sampleId = 0; sampleId < sampleSize; sampleId++) {
        currentSamples[sampleId] = combination(0, ranks.size() - 1, pathwaySize, gen);
        sort(currentSamples[sampleId].begin(), currentSamples[sampleId].end());

        double currentES = calcES(ranks, currentSamples[sampleId]);
        while (currentES <= 0) {
            vector<int> rnds = combination(0, ranks.size() - 1, pathwaySize, gen);
            sort(rnds.begin(), rnds.end());
            currentES = calcES(ranks, rnds);
            // rnds not saved to keep currentSamples uniform for ES+
            totalCount++;
        }
        posCount++;
        totalCount++;
    }

    posUnifScoreCount = make_pair(posCount, totalCount);

    duplicateSamples();
    while (ES > enrichmentScores.back()){
        for (int moves = 0; moves < sampleSize * pathwaySize;) {
            for (int sample_id = 0; sample_id < sampleSize; sample_id++) {
                moves += perturbate(ranks, currentSamples[sample_id], enrichmentScores.back(), gen);
            }
        }
        duplicateSamples();
        if (eps != 0){
            unsigned long k = enrichmentScores.size() / ((sampleSize + 1) / 2);
            if (k > - log2(0.5 * eps * exp(betaMeanLog(posUnifScoreCount.first, posUnifScoreCount.second)))) {
                break;
            }
        }
    }
}


pair<double, bool> EsRuler::getPvalue(double ES, double eps, bool sign) {
    unsigned long halfSize = (sampleSize + 1) / 2;

    auto it = enrichmentScores.begin();
    if (ES >= enrichmentScores.back()){
        it = enrichmentScores.end() - 1;
    }
    else{
        it = lower_bound(enrichmentScores.begin(), enrichmentScores.end(), ES);
    }

    unsigned long indx = 0;
    (it - enrichmentScores.begin()) > 0 ? (indx = (it - enrichmentScores.begin())) : indx = 0;

    unsigned long k = (indx) / halfSize;
    unsigned long remainder = sampleSize -  (indx % halfSize);

    double adjLog = betaMeanLog(halfSize, sampleSize);
    double adjLogPval = k * adjLog + betaMeanLog(remainder + 1, sampleSize);

    if (sign) {
        return make_pair(max(0.0, min(1.0, exp(adjLogPval))), true);
    } else {
        pair<double, bool> correction = calcLogCorrection(probCorrector, indx, posUnifScoreCount, sampleSize);
        double resLog = adjLogPval + correction.first;
        return make_pair(max(0.0, min(1.0, exp(resLog))), correction.second);
    }
}


int perturbate(const vector<double> &ranks, vector<int> &sample,
               double bound, mt19937 &rng) {
    double pertPrmtr = 0.1;
    int n = (int) ranks.size();
    int k = (int) sample.size();
    uniform_int_distribution<> uid_n(0, n - 1);
    uniform_int_distribution<> uid_k(0, k - 1);
    double NS = 0;
    for (int pos : sample) {
        NS += ranks[pos];
    }
    int iters = max(1, (int) (k * pertPrmtr));
    int moves = 0;
    for (int i = 0; i < iters; i++) {
        int id = uid_k(rng);
        int old = sample[id];
        NS -= ranks[sample[id]];

        sample[id] = uid_n(rng);
        while (id > 0 && sample[id] < sample[id - 1]) {
            swap(sample[id], sample[id - 1]);
            id--;
        }
        while (id < k - 1 && sample[id] > sample[id + 1]) {
            swap(sample[id], sample[id + 1]);
            id++;
        }

        if ((id > 0 && sample[id] == sample[id - 1]) || (id < k - 1 && sample[id] == sample[id + 1]) ||
            !compareStat(ranks, sample, NS + ranks[sample[id]], bound)) {
            // revert changes...
            sample[id] = old;
            while (id > 0 && sample[id] < sample[id - 1]) {
                swap(sample[id], sample[id - 1]);
                id--;
            }
            while (id < k - 1 && sample[id] > sample[id + 1]) {
                swap(sample[id], sample[id + 1]);
                id++;
            }
        } else {
            moves++;
        }
        NS += ranks[sample[id]];
    }
    return moves;
}
