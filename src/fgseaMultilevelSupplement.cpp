#include "fgseaMultilevelSupplement.h"
#include "esCalculation.h"

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

void fillRandomSample(set<int> &randomSample, mt19937 &gen,
                      const unsigned long ranksSize, const unsigned int pathwaySize) {
    randomSample.clear();
    uniform_int_distribution<> uid_n(0, static_cast<int>(ranksSize - 1));
    while (static_cast<int>(randomSample.size()) < pathwaySize) {
        randomSample.insert(uid_n(gen));
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

void EsRuler::extend(double ES, int seed, double eps) {
    unsigned int posCount = 0;
    unsigned int totalCount = 0;
    mt19937 gen(seed);

    for (int sampleId = 0; sampleId < sampleSize; sampleId++) {
        set<int> randomSample;
        fillRandomSample(randomSample, gen, ranks.size(), pathwaySize);
        currentSamples[sampleId] = vector<int>(randomSample.begin(), randomSample.end());
        double currentES = calcES(ranks, currentSamples[sampleId]);
        while (currentES <= 0) {
            fillRandomSample(randomSample, gen, ranks.size(), pathwaySize);
            currentES = calcES(ranks, vector<int>(randomSample.begin(), randomSample.end()));
            totalCount++;
        }
        posCount++;
        totalCount++;
    }

    posUnifScoreCount = make_pair(posCount, totalCount);

    duplicateSamples();
    while (ES > enrichmentScores.back()){
        for (int moves = 0; moves < sampleSize * pathwaySize;) {
            for (int sample_id = 0; sample_id < sampleSize; sample_id++) {
                moves += perturbate(ranks, currentSamples[sample_id], enrichmentScores.back(), gen);
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
