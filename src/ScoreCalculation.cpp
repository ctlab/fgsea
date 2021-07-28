#include "ScoreCalculation.h"

double getVarianceEstimator(const std::vector<double> & elems){
    double varHat = 0.0;
    for (auto elem : elems){
        varHat += elem * elem;
    }
    varHat = varHat / (elems.size() - 1);
    return varHat;
}


std::vector<double> getProfile(const std::vector<std::vector<double> > & E,
                               const std::vector<int> & indexes){
    unsigned m = E[0].size(); // number or rows must be nonzero
    std::vector<double> profile(m);

    for (unsigned j = 0; j < m; j++){
        double value = 0.0;
        for (auto indx : indexes){
            value += E[indx][j];
        }
        profile[j] = value;
    }
    return profile;
}


double getScore(const std::vector<double> & profile){
    return(getVarianceEstimator(profile));
}



void adjustProfile(const std::vector<std::vector<double> > & E,
                                  const std::vector<double> & profile,
                                  std::vector<double> & newProfile,
                                  unsigned idNew, unsigned idOld){
    for (unsigned i = 0; i < newProfile.size(); ++i){
        newProfile[i] = profile[i] - E[idOld][i] + E[idNew][i];
    }
}
