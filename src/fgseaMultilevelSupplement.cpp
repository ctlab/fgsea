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

EsRuler::SampleChunks::SampleChunks(int chunksNumber) : chunkSum(chunksNumber), chunkSize(chunksNumber), chunks(chunksNumber),
    chunkHasConvexHull(chunksNumber), chunkConvexHull(chunksNumber), chunkConvexHullBestPoint(chunksNumber) {}

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
    chunksNumber = max(1, (int) (0.5 * sqrt(pathwaySize)));
    chunkLastElement = vector<int>(chunksNumber);
    chunkLastElement[chunksNumber - 1] = ranks.size();
    vector<int> tmp(sampleSize);
    vector<SampleChunks> samplesChunks(sampleSize, SampleChunks(chunksNumber));

    duplicateSamples();
    while (ES > enrichmentScores.back()){
        for (int i = 0, pos = 0; i < chunksNumber - 1; ++i) {
            pos += (pathwaySize + i) / chunksNumber;
            for (int j = 0; j < sampleSize; ++j) {
                tmp[j] = currentSamples[j][pos];
            }
            nth_element(tmp.begin(), tmp.begin() + sampleSize / 2, tmp.end());
            chunkLastElement[i] = tmp[sampleSize / 2];
        }

        for (int i = 0; i < sampleSize; ++i) {
            fill(samplesChunks[i].chunkSum.begin(), samplesChunks[i].chunkSum.end(), 0.0);
            fill(samplesChunks[i].chunkSize.begin(), samplesChunks[i].chunkSize.end(), 0);
            int cnt = 0;
            samplesChunks[i].chunks[cnt].clear();
            for (int pos : currentSamples[i]) {
                while (chunkLastElement[cnt] <= pos) {
                    ++cnt;
                    samplesChunks[i].chunks[cnt].clear();
                }
                samplesChunks[i].chunks[cnt].push_back(pos);
                samplesChunks[i].chunkSum[cnt] += ranks[pos];
                samplesChunks[i].chunkSize[cnt]++;
            }
            fill(samplesChunks[i].chunkHasConvexHull.begin(), samplesChunks[i].chunkHasConvexHull.end(), false);
        }

        for (int moves = 0; moves < sampleSize * pathwaySize;) {
            for (int sampleId = 0; sampleId < sampleSize; sampleId++) {
                moves += perturbate(ranks, pathwaySize, samplesChunks[sampleId], enrichmentScores.back(), gen);
            }
        }

        for (int i = 0; i < sampleSize; ++i) {
            currentSamples[i].clear();
            for (int j = 0; j < chunksNumber; ++j) {
                for (int pos : samplesChunks[i].chunks[j]) {
                    currentSamples[i].push_back(pos);
                }
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

#define szof(x) ((int) (x).size())

int EsRuler::perturbate(const vector<double> &ranks, int k, EsRuler::SampleChunks &sampleChunks,
               double bound, mt19937 &rng) {
    double pertPrmtr = 0.1;
    int n = (int) ranks.size();
    uniform_int_distribution<> uid_n(0, n - 1);
    uniform_int_distribution<> uid_k(0, k - 1);
    double NS = 0;
    for (int i = 0; i < szof(sampleChunks.chunks); ++i) {
        for (int pos : sampleChunks.chunks[i]) {
            NS += ranks[pos];
        }
    }
    double q1 = 1.0 / (n - k);
    int iters = max(1, (int) (k * pertPrmtr));
    int moves = 0;

    for (int i = 0; i < iters; i++) {
        int oldInd = uid_k(rng);
        int oldChunkInd = 0, oldIndInChunk = 0;
        int oldVal;
        {
            int tmp = oldInd;
            while (szof(sampleChunks.chunks[oldChunkInd]) <= tmp) {
                tmp -= szof(sampleChunks.chunks[oldChunkInd]);
                ++oldChunkInd;
            }
            oldIndInChunk = tmp;
            oldVal = sampleChunks.chunks[oldChunkInd][oldIndInChunk];
        }

        int newVal = uid_n(rng);

        int newChunkInd = upper_bound(chunkLastElement.begin(), chunkLastElement.end(), newVal) - chunkLastElement.begin();
        int newIndInChunk = lower_bound(sampleChunks.chunks[newChunkInd].begin(), sampleChunks.chunks[newChunkInd].end(), newVal) - sampleChunks.chunks[newChunkInd].begin();

        if (newIndInChunk < szof(sampleChunks.chunks[newChunkInd]) && sampleChunks.chunks[newChunkInd][newIndInChunk] == newVal) {
            if (newVal == oldVal) {
                ++moves;
            }
            continue;
        }

        sampleChunks.chunks[oldChunkInd].erase(sampleChunks.chunks[oldChunkInd].begin() + oldIndInChunk);
        sampleChunks.chunks[newChunkInd].insert(
            sampleChunks.chunks[newChunkInd].begin() + newIndInChunk - (oldChunkInd == newChunkInd && oldIndInChunk < newIndInChunk ? 1 : 0), 
            newVal);

        NS = NS - ranks[oldVal] + ranks[newVal];

        sampleChunks.chunkSum[oldChunkInd] -= ranks[oldVal];
        sampleChunks.chunkSize[oldChunkInd]--;

        sampleChunks.chunkSum[newChunkInd] += ranks[newVal];
        sampleChunks.chunkSize[newChunkInd]++;

        sampleChunks.chunkHasConvexHull[oldChunkInd] = false;
        sampleChunks.chunkHasConvexHull[newChunkInd] = false;

        double q2 = 1.0 / NS;

        double val = 0;
        int last = -1;
        int lastChunkLastElement = 0;

        bool ok = false;

        for (int i = 0; i < chunksNumber; ++i) {
            if (!szof(sampleChunks.chunks[i]) || val + q2 * sampleChunks.chunkSum[i] < bound) {
                val += q2 * sampleChunks.chunkSum[i] - q1 * (chunkLastElement[i] - lastChunkLastElement - sampleChunks.chunkSize[i]);
                last = chunkLastElement[i] - 1;
                lastChunkLastElement = chunkLastElement[i];
                continue;
            }
            auto& convexHull = sampleChunks.chunkConvexHull[i];
            int& posConvexHull = sampleChunks.chunkConvexHullBestPoint[i];
            if (sampleChunks.chunkHasConvexHull[i]) {
                double cval = -convexHull[posConvexHull].first * q1 + convexHull[posConvexHull].second * q2;
                while (0 < posConvexHull && val + cval <= bound) {
                    double nval = -convexHull[posConvexHull - 1].first * q1 + convexHull[posConvexHull - 1].second * q2;
                    if (nval > cval) {
                        cval = nval;
                        --posConvexHull;
                    } else {
                        break;
                    }
                }
                if (val + cval > bound) {
                    ok = true;
                    break;
                }
                while (posConvexHull < szof(convexHull) - 1 && val + cval <= bound) {
                    double nval = -convexHull[posConvexHull + 1].first * q1 + convexHull[posConvexHull + 1].second * q2;
                    if (nval > cval) {
                        cval = nval;
                        ++posConvexHull;
                    } else {
                        break;
                    }
                }
                if (val + cval > bound) {
                    ok = true;
                    break;
                }

                val += q2 * sampleChunks.chunkSum[i] - q1 * (chunkLastElement[i] - lastChunkLastElement - sampleChunks.chunkSize[i]);
                last = chunkLastElement[i] - 1;
            } else {
                int curX = 0;
                double curY = 0;
                convexHull.clear();
                convexHull.push_back({0, 0});
                double best = -1e100;
                for (int pos : sampleChunks.chunks[i]) {
                    curY += ranks[pos];
                    curX += pos - last - 1;
                    val += ranks[pos] * q2 - (pos - last - 1) * q1;
                    if (val > bound) {
                        ok = true;
						// break;
                    }
                    last = pos;
                    while (szof(convexHull) >= 2) {
                        int dx1 = convexHull.back().first - convexHull[szof(convexHull) - 2].first;
                        double dy1 = convexHull.back().second - convexHull[szof(convexHull) - 2].second;
                        int dx2 = curX - convexHull.back().first;
                        double dy2 = curY - convexHull.back().second;
                        if (dx1 * dy2 - dx2 * dy1 >= 0) {
                            convexHull.pop_back();
                        } else {
                            break;
                        }
                    }
                    convexHull.push_back({curX, curY});
                    if (val > best) {
                        best = val;
                        posConvexHull = szof(convexHull) - 1;
                    }
                }
                curX += chunkLastElement[i] - last - 1;
                while (szof(convexHull) >= 2) {
                    int dx1 = convexHull.back().first - convexHull[szof(convexHull) - 2].first;
                    double dy1 = convexHull.back().second - convexHull[szof(convexHull) - 2].second;
                    int dx2 = curX - convexHull.back().first;
                    double dy2 = curY - convexHull.back().second;
                    if (dx1 * dy2 - dx2 * dy1 >= 0) {
                        convexHull.pop_back();
                    } else {
                        break;
                    }
                }
                convexHull.push_back({curX, curY});
                sampleChunks.chunkHasConvexHull[i] = true;
                val -= q1 * (chunkLastElement[i] - last - 1);
                last = chunkLastElement[i] - 1;
                if (ok) {
                    break;
                }
            }
            lastChunkLastElement = chunkLastElement[i];
        }

        if (!ok) {
        	NS = NS - ranks[newVal] + ranks[oldVal];
            
            sampleChunks.chunkSum[oldChunkInd] += ranks[oldVal];
            sampleChunks.chunkSize[oldChunkInd]++;

            sampleChunks.chunkSum[newChunkInd] -= ranks[newVal];
            sampleChunks.chunkSize[newChunkInd]--;

            sampleChunks.chunks[newChunkInd].erase(
                sampleChunks.chunks[newChunkInd].begin() + newIndInChunk - (oldChunkInd == newChunkInd && oldIndInChunk < newIndInChunk ? 1 : 0));
            sampleChunks.chunks[oldChunkInd].insert(sampleChunks.chunks[oldChunkInd].begin() + oldIndInChunk, oldVal);

            sampleChunks.chunkHasConvexHull[oldChunkInd] = false;
            sampleChunks.chunkHasConvexHull[newChunkInd] = false;
        } else {
            ++moves;
        }
    }
    return moves;
}