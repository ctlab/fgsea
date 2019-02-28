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


double cross(int x1, double y1, int x2, double y2) {
  return x1 * y2 - y1 * x2;
}

struct genes {
  vector<double> const& ranks;
  int chunk_length;
  vector<int> inds;
  vector<int> next;
  vector<int> first_in_chunk, size_of_chunk;
  vector<double> sum_in_chunk;
  int diagonal_x;
  double diagonal_y;

  genes(vector<double> const& _ranks, vector<int> const& _inds) : ranks(_ranks), chunk_length((int) (ranks.size() / sqrt(_inds.size()))), 
    inds(_inds), next(inds.size(), -1), first_in_chunk((ranks.size() + chunk_length - 1) / chunk_length, -1), size_of_chunk(first_in_chunk.size()), sum_in_chunk(first_in_chunk.size()),
    diagonal_x(ranks.size() - inds.size()), diagonal_y(0) {
    int chunk = -1;
    for (int i = 0; i < inds.size(); ++i) {
      diagonal_y += ranks[inds[i]];
      int cur = inds[i] / chunk_length;
      size_of_chunk[cur]++;
      sum_in_chunk[cur] += ranks[inds[i]];
      if (cur > chunk) {
        chunk = cur;
        first_in_chunk[chunk] = i;
      } else {
        next[i - 1] = i;
      }
    }
  }

  // void remove_pos(int, int);
  // void add_ind(int, int, int);
  // int add_ind(int, int);
  bool check_es_not_less(double bound) {
    bound *= diagonal_x * diagonal_y;
    double y = 0;
    int x = 0;
    int prev = 0;
    for (int i = 0; i < first_in_chunk.size(); ++i) {
      if (cross(diagonal_x, diagonal_y, x, y + sum_in_chunk[i]) <= bound) {
        x += chunk_length * (i + 1) - prev - size_of_chunk[i];
        y += sum_in_chunk[i];
        prev = (i + 1) * chunk_length;
        continue;
      }
      int pos = first_in_chunk[i];
      while (pos != -1) {
        x += inds[pos] - prev;
        y += ranks[inds[pos]];
        prev = inds[pos] + 1;
        pos = next[pos];
        if (cross(diagonal_x, diagonal_y, x, y) > bound) {
          return true;
        }
      }
    }

    return false;
  }

  vector<int> get_inds() {
    vector<int> ret;
    for (int i = 0; i < first_in_chunk.size(); ++i) {
      int pos = first_in_chunk[i];
      while (pos != -1) {
        ret.push_back(inds[pos]);
        pos = next[pos];
      }
    }
    return ret;
  }

  bool try_replace(int order_number, int new_ind, double bound) {
    int pos = -1;
    int old_prev = -1;
    for (int i = 0; i < first_in_chunk.size(); ++i) {
      if (order_number < size_of_chunk[i]) {
        pos = first_in_chunk[i];
        while (order_number) {
          old_prev = pos;
          pos = next[pos];
          --order_number;
        }
        break;
      }
      order_number -= size_of_chunk[i];
    }
    int old_ind = inds[pos];
    if (old_ind == new_ind) {
      return false;
    }

    diagonal_y += -ranks[old_ind] + ranks[new_ind];

    int old_chunk = old_ind / chunk_length;
    int new_chunk = new_ind / chunk_length;

    // remove_pos(old_prev, pos);
    size_of_chunk[old_chunk]--;
    sum_in_chunk[old_chunk] -= ranks[old_ind];
    if (old_prev == -1) {
      first_in_chunk[old_chunk] = next[pos];
    } else {
      next[old_prev] = next[pos];
    }


    // int new_prev = add_ind(new_ind, pos);
    inds[pos] = new_ind;
    size_of_chunk[new_chunk]++;
    sum_in_chunk[new_chunk] += ranks[new_ind];
    int new_prev = -1;
    int tmp = first_in_chunk[new_chunk];
    while (tmp != -1 && inds[tmp] <= new_ind) {
      new_prev = tmp;
      tmp = next[tmp];
    }

    if (new_prev == -1) {
      next[pos] = first_in_chunk[new_chunk];
      first_in_chunk[new_chunk] = pos;
    } else {
      next[pos] = next[new_prev];
      next[new_prev] = pos;
    }

    if ((new_prev == -1 || inds[new_prev] < new_ind) && check_es_not_less(bound)) {
      return true;
    }

    // remove_pos(new_prev, pos);
    size_of_chunk[new_chunk]--;
    sum_in_chunk[new_chunk] -= ranks[new_ind];
    if (new_prev == -1) {
      first_in_chunk[new_chunk] = next[pos];
    } else {
      next[new_prev] = next[pos];
    }


    // add_ind(old_ind, pos, old_prev);
    inds[pos] = old_ind;
    size_of_chunk[old_chunk]++;
    sum_in_chunk[old_chunk] += ranks[old_ind];

    if (old_prev == -1) {
      next[pos] = first_in_chunk[old_chunk];
      first_in_chunk[old_chunk] = pos;
    } else {
      next[pos] = next[old_prev];
      next[old_prev] = pos;
    }

    diagonal_y += ranks[old_ind] - ranks[new_ind];

    return false;
  }
};

int perturbate(const vector<double> &S, vector<int> &p, double bound, mt19937& rng, double pert_coeff) {
  int n = (int) S.size();
  int k = (int) p.size();
  genes gns(S, p);
  uniform_int_distribution<> uid_n(0, n - 1);
  uniform_int_distribution<> uid_k(0, k - 1);
  double NS = 0;
  for (int pos : p) {
    NS += S[pos];
  }
  int moves = 0;
  int iters = max(1, (int) (k * pert_coeff));
  for (int i = 0; i < iters; i++) {
    int pos = uid_k(rng);
    int new_ind = uid_n(rng);
    moves += gns.try_replace(pos, new_ind, bound);
  }

  p = gns.get_inds();
  return moves;
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


  for (int sample_id = 0; 2 * sample_id < sampleSize; sample_id++) {
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
      for (int moves = 0; moves < sampleSize * pathwaySize; ){
         for (int sample_id = 0; sample_id < sampleSize; sample_id++) {
           moves += perturbate(S, esPvalObj.sets[sample_id], esPvalObj.cutoffs.back(), gen, 0.1);
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
    double probStatPos = exp(boost::math::digamma(esPvalObj.posStatNum) - boost::math::digamma(sampleSize + 1));

    auto it = lower_bound(esPvalObj.cutoffs.begin(), esPvalObj.cutoffs.end(), enrichmentScore);
    int k = ((it - esPvalObj.cutoffs.begin())) / halfSize;
    int remainder = sampleSize - (it - esPvalObj.cutoffs.begin()) % (halfSize);
    double adjLog = boost::math::digamma(halfSize) - boost::math::digamma(sampleSize + 1);
    double adjLogPval = k * adjLog + (boost::math::digamma(remainder) - boost::math::digamma(sampleSize + 1));
    pval = exp(adjLogPval);

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
