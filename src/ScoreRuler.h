#ifndef SCORE_RULER_H
#define SCORE_RULER_H

#include <vector>
#include <random>
#include <algorithm>
#include "ScoreCalculation.h"



class ScoreRuler {
private:
    const std::vector<std::vector<double> > & expressionMatrix;
    const unsigned sampleSize;
    const unsigned genesetSize;

    std::vector<double> scores;
    std::vector<std::vector<unsigned> > currentSample;
    std::vector<std::vector<double> > currentProfiles;


    void duplicateSampleElements();
    int updateElement(std::vector<unsigned> & element, std::vector<double> & profile,
                      double threshold, std::mt19937 &mtGen);
public:
    ScoreRuler(const std::vector<std::vector<double> > & inpE,
               unsigned inpSampleSize, unsigned inpGenesetSize);
    ~ScoreRuler();

    void extend(double inpScore, int seed, double eps);
    double getPvalue(double inpScore, double eps);

    // pair<double, bool> getPvalue(double ES, double eps, bool sign);
};

#endif //SCORE_RULER_H
