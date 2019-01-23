#include "fgseaMultilevelSupplement.h"

#include <random>
#include <vector>
#include <cmath>
#include <algorithm>
#include <set>
using namespace std;


#include <Rcpp.h>
using namespace Rcpp;

double calcES(const vector<double>& S, const vector<int>& p, double NS) {
  // p must be sorted
  int n = (int) S.size();
  int k = (int) p.size();
  double res = 0.0;
  double cur = 0.0;
  double q1 = 1.0 / (n - k);
  double q2 = 1.0 / NS;
  int last = -1;
  for (int pos : p) {
    cur -= q1 * (pos - last - 1);
    if (abs(cur) > abs(res)) {
      res = cur;
    }
    cur += q2 * S[pos];
    if (abs(cur) > abs(res)) {
      res = cur;
    }
    last = pos;
  }
  return res;
}

double calcES(const vector<double>& S, const vector<int>& p) {
  // p must be sorted
  double NS = 0.0;
  for (int pos : p) {
    NS += S[pos];
  }
  return calcES(S, p, NS);
}

double calcPositiveES(const vector<double>& S, const vector<int>&p, double NS) {
  // p must be sorted
  int n = (int) S.size();
  int k = (int) p.size();
  double res = 0.0;
  double cur = 0.0;
  double q1 = 1.0 / (n - k);
  double q2 = 1.0 / NS;
  int last = -1;
  for (int pos : p) {
    cur += q2 * S[pos] - q1 * (pos - last - 1);
    res = max(res, cur);
    last = pos;
  }
  return res;
}

double calcPositiveES(const vector<double>& S, const vector<int>& p) {
  // p must be sorted
  double NS = 0.0;
  for (int pos : p) {
    NS += S[pos];
  }
  return calcPositiveES(S, p, NS);
}

int compareStat(const vector<double>& S, const vector<int>& p, double NS, double bound) {
  // p must be sorted
  int n = (int) S.size();
  int k = (int) p.size();
  double cur = 0.0;
  double q1 = 1.0 / (n - k);
  double q2 = 1.0 / NS;
  int last = -1;
  for (int pos : p) {
    cur += q2 * S[pos] - q1 * (pos - last - 1);
    if (cur > bound) {
      return 1;
    }
    last = pos;
  }
  return -1;
}

void perturbate(const vector<double> &S, vector<int> &p, double bound, mt19937& rng, double pert_coeff) {
  int n = (int) S.size();
  int k = (int) p.size();
  uniform_int_distribution<> uid_n(0, n - 1);
  uniform_int_distribution<> uid_k(0, k - 1);
  double NS = 0;
  for (int pos : p) {
    NS += S[pos];
  }
  int iters = max(1, (int) (k * pert_coeff));
  for (int i = 0; i < iters; i++) {
    int id = uid_k(rng);
    int old = p[id];
    NS -= S[p[id]];
    p[id] = uid_n(rng);
    while (id > 0 && p[id] < p[id - 1]) {
      swap(p[id], p[id - 1]);
      id--;
    }
    while (id < k - 1 && p[id] > p[id + 1]) {
      swap(p[id], p[id + 1]);
      id++;
    }
    if ((id > 0 && p[id] == p[id - 1]) || (id < k - 1 && p[id] == p[id + 1]) || compareStat(S, p, NS + S[p[id]], bound) == -1) {
      // revert changes...
      p[id] = old;
      while (id > 0 && p[id] < p[id - 1]) {
        swap(p[id], p[id - 1]);
        id--;
      }
      while (id < k - 1 && p[id] > p[id + 1]) {
        swap(p[id], p[id + 1]);
        id++;
      }
    }
    NS += S[p[id]];
  }
}

void sort_by_scores(const vector<double>& S, vector< vector<int> >& sets)
{
  vector< pair<double, int>> scores_pos;
  for (int i = 0; i < sets.size(); i++)
  {
    scores_pos.push_back(make_pair(calcPositiveES(S, sets[i]), i));
  }
  sort(scores_pos.begin(), scores_pos.end());

  vector<vector<int> > res(sets.size());

  for (int i = 0; i < sets.size(); i++) {
      swap(res[i], sets[scores_pos[i].second]);
   }
  swap(res, sets);
}

bool check_medians(const vector<double>& S,const vector<vector<int> >& sets)
{
  vector<double> ES_vector(sets.size());
  for (int i = 0; i < ES_vector.size(); i++)
  {
    ES_vector[i] = calcPositiveES(S, sets[i]);
  }
  if (ES_vector.size() > 2)
  {
    auto middle_it = ES_vector.begin() + (ES_vector.size()/2);
    vector<double> left_part(ES_vector.begin(), middle_it - 1);
    vector<double> right_part(middle_it + 1, ES_vector.end());

    nth_element(left_part.begin(),
                left_part.begin() + left_part.size() / 2,
                left_part.end());
     nth_element(right_part.begin(),
                 right_part.begin() + right_part.size() / 2,
                 right_part.end());

    double left_median = left_part[left_part.size() / 2];
    double right_median = right_part[right_part.size() / 2]; 
    return (left_median > right_median);
  }
  else if (ES_vector.size() == 2)
  {
    return ES_vector[0] > ES_vector[1];
  }
  else return false;
}

void duplicateSets(EsPvalConnection &esPvalObj, int sampleSize, const vector<double> &S)
{
  // Duplicate sets of genes with enrichment score greater than the median value.
  int posStatCount = 0;
  vector<pair<double, int> > stats(sampleSize);
  for (int sample_id = 0; sample_id < sampleSize; sample_id++)
  {
      double stat = calcPositiveES(S, esPvalObj.sets[sample_id]);
      double stat_real = calcES(S, esPvalObj.sets[sample_id]);
      if (esPvalObj.cutoffs.empty())
      {
          esPvalObj.random_pairs.emplace_back(stat, stat_real);
      }
      if (stat_real > 0)
      {
          posStatCount++;
      }
      stats[sample_id] = make_pair(stat, sample_id);
  }
  sort(stats.begin(), stats.end());
  if (esPvalObj.posStatNum == 0)
  {
      esPvalObj.posStatNum = posStatCount;
  }


  for (int sample_id = 0; 2 * sample_id < sampleSize - 1; sample_id++) {
    esPvalObj.cutoffs.emplace_back(stats[sample_id].first);
  }



  vector< vector<int> > new_sets;
  for (int sample_id = 0; 2 * sample_id < sampleSize - 2; sample_id++) {
      for (int rep = 0; rep < 2; rep++) {
          new_sets.push_back(esPvalObj.sets[stats[sampleSize - 1 - sample_id].second]);
      }
  }

  new_sets.push_back(esPvalObj.sets[stats[sampleSize >> 1].second]);
  swap(esPvalObj.sets, new_sets);
  sort_by_scores(S, esPvalObj.sets);
  return;
}

void calcPvalues(EsPvalConnection &esPvalObj, vector<double> S, int pathwaySize,
                 double ES, int sampleSize, int seed, double absEps)
{
    mt19937 gen(seed);
    uniform_int_distribution<> uid_n(0, S.size() - 1);
    for (int sample_id = 0; sample_id < sampleSize; sample_id++)
    {

        set<int> random_set;
        while ((int) random_set.size() < pathwaySize)
        {
            random_set.insert(uid_n(gen));
        }
        esPvalObj.sets[sample_id] = vector<int>(random_set.begin(), random_set.end());
    }

    duplicateSets(esPvalObj, sampleSize, S);

    while (true)
    {
      int k = 2 * (esPvalObj.cutoffs.size() / (sampleSize + 1));
      if (ES < esPvalObj.cutoffs.back() || k > -log2(absEps)){
        break;
      }
      bool median_checker = false;
      int  todo_num = 0;
      for (int i = 0; (!median_checker) || (i < todo_num); i++){
        if (check_medians(S, esPvalObj.sets) && todo_num == 0){
          median_checker = true;
          todo_num = i*4;
        }
        for (int sample_id = 0; sample_id < sampleSize; sample_id++) {
            perturbate(S, esPvalObj.sets[sample_id], esPvalObj.cutoffs.back(), gen, 0.1);
        }
      }
      duplicateSets(esPvalObj, sampleSize, S);
    }
    return;
}


double findEsPval(const EsPvalConnection &esPvalObj, double enrichmentScore, int sampleSize, bool sign)
{
    int halfSize = (sampleSize + 1) / 2;
    double pval = 0;
    double probStatPos = pow(M_E, boost::math::digamma(esPvalObj.posStatNum) - boost::math::digamma(sampleSize + 1));
    if (enrichmentScore < esPvalObj.cutoffs.front())
    {
        pval = pow(M_E, boost::math::digamma(sampleSize) - boost::math::digamma(sampleSize + 1));
    }
    else if (enrichmentScore > esPvalObj.cutoffs.back())
    {
      int k = (esPvalObj.cutoffs.size() / halfSize);
      int remainder = sampleSize - (esPvalObj.cutoffs.size()) % (halfSize);
      double adjLog = boost::math::digamma(halfSize) - boost::math::digamma(sampleSize + 1);
      double adjLogPval = k * adjLog + (boost::math::digamma(remainder) - boost::math::digamma(sampleSize + 1));
      pval = pow(M_E, adjLogPval);
    }
    else
    {
      auto it = lower_bound(esPvalObj.cutoffs.begin(), esPvalObj.cutoffs.end(), enrichmentScore);
      int k = ((it - esPvalObj.cutoffs.begin())) / halfSize;
      int remainder = sampleSize - (it - esPvalObj.cutoffs.begin()) % (halfSize);
      double adjLog = boost::math::digamma(halfSize) - boost::math::digamma(sampleSize + 1);
      double adjLogPval = k * adjLog + (boost::math::digamma(remainder) - boost::math::digamma(sampleSize + 1));
      pval = pow(M_E, adjLogPval);
    }
    if (sign){
      pval = max(0.0, min(1.0, pval));
    }
    else{
      double correction;
      int badSets = 0;
      int totalSets = 0;
      for (auto &pp : esPvalObj.random_pairs)
      {
          if (pp.second <= enrichmentScore)
          {
              badSets += (pp.first > enrichmentScore);
          }
          totalSets++;
      }
      correction = badSets*1.0 / totalSets;
      pval = max(0.0, min(1.0, (pval - correction)/probStatPos));
    }
    return pval;
}
