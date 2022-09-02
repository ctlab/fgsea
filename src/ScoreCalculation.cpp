#include "ScoreCalculation.h"

double getVarianceEstimator(const std::vector<double> & elems){
    double varHat = 0.0;
    for (double elem : elems){
        varHat += elem * elem;
    }
    varHat = varHat / (elems.size() - 1);
    return varHat;
}


std::vector<double> getProfile(const std::vector<double> & E,
                               const std::vector<unsigned> & indexes,
                               unsigned m){
    // unsigned m = E[0].size(); // number or rows must be nonzero
    std::vector<double> profile(m);

    for (unsigned j = 0; j < m; j++){
        double value = 0.0;
        for (unsigned indx : indexes){
            value += E[indx*m + j];
        }
        profile[j] = value;
    }
    return profile;
}


double getScore(const std::vector<double> & profile){
    return(getVarianceEstimator(profile));
}



void adjustProfile(const std::vector<double> & E,
                                  const std::vector<double> & profile,
                                  std::vector<double> & newProfile,
                                  unsigned idNew, unsigned idOld,
                                  unsigned m){
    for (unsigned i = 0; i < newProfile.size(); i++){
        newProfile[i] = profile[i] - E[idOld*m + i] + E[idNew*m + i];
    }
}
