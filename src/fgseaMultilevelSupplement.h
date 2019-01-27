#pragma once

#include <random>
#include <vector>
using namespace std;

// [[Rcpp::depends(BH)]]
#include <Rcpp.h>
using namespace Rcpp;

#include <boost/math/special_functions/digamma.hpp>

struct EsPvalConnection
{
    vector<vector<int> > sets;
    vector<pair<double, double> > random_pairs;
    vector<double> cutoffs;
    int posStatNum;
    EsPvalConnection(int setsSize)
    {
        this -> sets.resize(setsSize);
        this -> posStatNum = 0;
    }
};

// [[Rcpp::plugins(cpp17)]]
void calcPvalues(EsPvalConnection &esPvalObj, vector<double> S, int pathwaySize,
                 double ES, int sampleSize, int seed, double absEps);

// [[Rcpp::plugins(cpp17)]]
double findEsPval(const EsPvalConnection &esPvalObj, double enrichmentScore, int sampleSize, bool sign);

// [[Rcpp::plugins(cpp17)]]
double correctionPval(double ES, double pval, vector<pair<double, double> > randomPairs);


double calcES(const vector<double>& S, const vector<int>& p, const double NS);
double calcES(const vector<double>& S, const vector<int>& p);

double calcPositiveES(const vector<double>& S, const vector<int>&p, const double NS);
double calcPositiveES(const vector<double>& S, const vector<int>& p);

int compareStat(const vector<double>& S, const vector<int>& p, const double NS, const double bound);

int perturbate(const vector<double> &S, vector<int> &p, double bound, mt19937& rng, double pert_coeff);