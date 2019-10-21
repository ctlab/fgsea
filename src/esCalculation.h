#ifndef FGSEAMULTILEVELCPP_ESCALCULATION_H
#define FGSEAMULTILEVELCPP_ESCALCULATION_H

#include <vector>
#include <cmath>

using namespace std;

double calcES(const vector<double> &ranks, const vector<int> &p, double NS);

double calcES(const vector<double> &ranks, const vector<int> &p);

double calcPositiveES(const vector<double> &ranks, const vector<int> &p, double NS);

double calcPositiveES(const vector<double> &ranks, const vector<int> &p);

bool compareStat(const vector<double> &ranks, const vector<int> &p, double NS, double bound);

#endif //FGSEAMULTILEVELCPP_ESCALCULATION_H
