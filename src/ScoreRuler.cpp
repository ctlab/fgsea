#include "ScoreRuler.h"
#include "util.h"



ScoreRuler::ScoreRuler(const std::vector<std::vector<double> > & inpE,
                       unsigned inpSampleSize, unsigned inpGenesetSize):
    expressionMatrix(inpE), sampleSize(inpSampleSize), genesetSize(inpGenesetSize){
    currentSample.resize(inpSampleSize), currentProfiles.resize(inpSampleSize);
}



ScoreRuler::~ScoreRuler() = default;


void ScoreRuler::duplicateSampleElements(){
    /*
     * Replaces sample elements that are less than median with elements
     * that are greater that the median
     */

    std::vector<std::pair<double, unsigned> > scoreAndIndex(sampleSize);

    for (unsigned elemIndex = 0; elemIndex < sampleSize; elemIndex++) {
        double elemScore = getScore(currentProfiles[elemIndex]);
        scoreAndIndex[elemIndex] = std::make_pair(elemScore, elemIndex);
    }
    std::sort(scoreAndIndex.begin(), scoreAndIndex.end());

    for (unsigned elemIndex = 0; 2 * elemIndex < sampleSize; elemIndex++){
        scores.push_back(scoreAndIndex[elemIndex].first);
    }

    std::vector<std::vector<int> > tempSample;
    std::vector<std::vector<double> > tempProfiles;
    for (unsigned elemIndex = 0; 2 * elemIndex < sampleSize - 2; elemIndex++) {
        for (unsigned rep = 0; rep < 2; rep++) {
            tempSample.push_back(currentSample[scoreAndIndex[sampleSize - 1 - elemIndex].second]);
            tempProfiles.push_back(currentProfiles[scoreAndIndex[sampleSize - 1 - elemIndex].second]);
        }
    }
    tempSample.push_back(currentSample[scoreAndIndex[sampleSize >> 1].second]);
    tempProfiles.push_back(currentProfiles[scoreAndIndex[sampleSize >> 1].second]);
    std::swap(currentSample, tempSample);
    std::swap(currentProfiles, tempProfiles);
}

void ScoreRuler::extend(double inpScore, int seed, double eps) {
    std::mt19937 mtGen(seed);

    std::vector<unsigned> indexSpace(expressionMatrix.size());
    std::iota(indexSpace.begin(), indexSpace.end(), 0);
    // fill currentSample
    for (int elemIndex = 0; elemIndex < sampleSize; elemIndex++) {
        currentSample[elemIndex] = combination(0, indexSpace.size() - 1, genesetSize, mtGen);
        currentProfiles[elemIndex] = getProfile(expressionMatrix, currentSample[elemIndex]);
    }

    duplicateSampleElements();
    while (scores.back() <= inpScore - 1e-10){
        for (int moves = 0; moves < sampleSize * genesetSize;) {
            for (int elemIndex = 0; elemIndex < sampleSize; elemIndex++) {
                moves += updateElement(currentSample[elemIndex], currentProfiles[elemIndex], scores.back(), mtGen);
            }
        }

        duplicateSampleElements();
        if (eps != 0){
            unsigned long k = scores.size() / ((sampleSize + 1) / 2);
            if (k > - log2(0.5 * eps)) {
                break;
            }
        }
    }
}

double ScoreRuler::getPvalue(double inpScore, double eps){
    unsigned long halfSize = (sampleSize + 1) / 2;

    auto it = scores.begin();
    if (inpScore >= scores.back()){
        it = scores.end() - 1;
    }
    else{
        it = lower_bound(scores.begin(), scores.end(), inpScore);
    }

    unsigned long indx = 0;
    (it - scores.begin()) > 0 ? (indx = (it - scores.begin())) : indx = 0;

    unsigned long k = (indx) / halfSize;
    unsigned long remainder = sampleSize -  (indx % halfSize);

    double adjLog = betaMeanLog(halfSize, sampleSize);
    double adjLogPval = k * adjLog + betaMeanLog(remainder + 1, sampleSize);
    return std::max(0.0, std::min(1.0, exp(adjLogPval)));
}


int ScoreRuler::updateElement(std::vector<int> & element,
                              std::vector<double> & profile,
                              double threshold,
                              std::mt19937 &mtGen){
    double upPrmtr = 0.1;
    auto n = (unsigned) (expressionMatrix.size());

    uid_wrapper uid_n(0, n - 1, mtGen);
    uid_wrapper uid_k(0, genesetSize - 1, mtGen);

    unsigned niters = std::max(1, (int) (genesetSize * upPrmtr));
    int moves = 0;
    std::vector<double> newProfile(profile.size());
    for (unsigned i = 0; i < niters; i++){
        unsigned toDrop = uid_k();
        unsigned indxOld = element[toDrop];
        unsigned indxNew = uid_n();

        adjustProfile(expressionMatrix, profile, newProfile, indxNew, indxOld);
        double newScore = getScore(newProfile);

        if (newScore < threshold || (*find(element.begin(), element.end(), indxNew) == indxNew)){
            element[toDrop] = indxOld;
        } else{
            element[toDrop] = indxNew;
            profile.swap(newProfile);
            moves++;
        }
    }
    return moves;
}
