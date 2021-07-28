#ifndef SCORE_CALC_H
#define SCORE_CALC_H
#include <cmath>
#include <algorithm>
#include <numeric>
#include <vector>

double getVarianceEstimator(const std::vector<double> & elems);

double getScore(const std::vector<double> & profile);

std::vector<double> getProfile(const std::vector<std::vector<double> > & E, const std::vector<int> & indexes);

void adjustProfile(const std::vector<std::vector<double> > & E,
                  const std::vector<double> & profile,
                  std::vector<double> & newProfile,
                  unsigned idNew, unsigned idOld);

#endif //SCORE_CALC_H
