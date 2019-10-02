#include <iostream>
#include <string>
#include <algorithm>
#include <set>
#include <vector>
#include <map>
#include <fstream>
using namespace std;


int main() {
  freopen("roundRanks.txt", "r", stdin);
  int n;
  scanf("%d", &n);
  vector <int> rank_rounded(n);
  for (int i = 0; i < n; i++) {
    scanf("%d", &rank_rounded[i]);
  }
  fclose(stdin);
  freopen("inpPathways.txt", "r", stdin);
  int test_count;
  scanf("%d", &test_count);
  map<pair<int,double>,pair<double,double>> memo;

  ofstream out_f;
  out_f.open("exactResults.tsv");
  out_f << "size\tES\tp_exact\ttotal_removed" << endl;
  for (int test_id = 0; test_id < test_count; test_id++) {
    int k;
    double gsea;
    scanf("%d %lf", &k, &gsea);
    cout << "Current request: " << test_id + 1 << "/" << test_count;
    cout << ". (Pathway Size: " << k << ", ES: " << gsea << ")" << endl;
    double ans = 0.0;
    double eps = 1e-40;
    double total_removed = 0.0;
    if (memo.find(make_pair(k, gsea)) != memo.end()) {
      ans = memo[make_pair(k, gsea)].first;
      total_removed = memo[make_pair(k, gsea)].second;
    } else {
      int flag = 0;
      if (gsea < 0) {
        reverse(rank_rounded.begin(), rank_rounded.end());
        gsea = -gsea;
        flag = 1;
      }
      if (k <= 500) {
        vector < vector < vector <double> > > dp(k + 1);
        vector < vector <int> > dp_from(k + 1);
        vector < vector <double> > easy(k + 1);
        dp[0].push_back(vector <double>(1, 1.0));
        dp_from[0].push_back(0);
        int milestone = 0;
        for (int i = 0; i < n; i++) {
          if (i >= milestone) {
            milestone += n / 10;
          }
          for (int j = k - 1; j >= 0; j--) {
            if ((i - j) * 1.0 / (n - k) + gsea > 1.0) {
              for (int s = 0; s < (int) dp[j].size(); s++) {
                for (int bound = dp_from[j][s]; bound < dp_from[j][s] + (int) dp[j][s].size(); bound++) {
                  if ((int) easy[j].size() <= bound - s + 1) {
                    easy[j].resize(bound - s + 2);
                  }
                  easy[j][bound - s + 1] += dp[j][s][bound - dp_from[j][s]];
                }
              }
              dp_from[j].clear();
              dp[j].clear();
            }
            double p_take = (k - j) * 1.0 / (n - i);
            for (int s = 0; s < (int) dp[j].size(); s++) {
              while (!dp[j][s].empty() && dp[j][s].back() < eps) {
                total_removed += dp[j][s].back();
                dp[j][s].pop_back();
              }
              while (!dp[j][s].empty() && dp[j][s][0] < eps) {
                total_removed += dp[j][s][0];
                dp[j][s].erase(dp[j][s].begin());
                dp_from[j][s]++;
              }
              if (dp[j][s].empty()) {
                continue;
              }
              int new_s = s + rank_rounded[i];
              while ((int) dp[j + 1].size() <= new_s) {
                dp[j + 1].push_back(vector <double>());
                dp_from[j + 1].push_back(1 << 30);
              }
              int aux_bound = max(new_s - 1, (int) (new_s / ((i - j) * 1.0 / (n - k) + gsea)));
              int min_new_bound = max(dp_from[j][s], aux_bound);
              if (min_new_bound < dp_from[j + 1][new_s]) {
                if (!dp[j + 1][new_s].empty()) {
                  dp[j + 1][new_s].insert(dp[j + 1][new_s].begin(), dp_from[j + 1][new_s] - min_new_bound, 0.0);
                }
                dp_from[j + 1][new_s] = min_new_bound;
              }
              int max_new_bound = max(dp_from[j][s] + (int) dp[j][s].size() - 1, aux_bound);
              if (max_new_bound >= dp_from[j + 1][new_s] + (int) dp[j + 1][new_s].size()) {
                dp[j + 1][new_s].resize(max_new_bound - dp_from[j + 1][new_s] + 1);
              }
              for (int bound = dp_from[j][s]; bound < dp_from[j][s] + (int) dp[j][s].size(); bound++) {
                int new_bound = max(bound, aux_bound);
                dp[j + 1][new_s][new_bound - dp_from[j + 1][new_s]] += dp[j][s][bound - dp_from[j][s]] * p_take;
                dp[j][s][bound - dp_from[j][s]] *= 1.0 - p_take;
              }
            }
            for (int s = 0; s < (int) easy[j].size(); s++) {
              int new_s = max(0, s - rank_rounded[i]);
              if ((int) easy[j + 1].size() <= new_s) {
                easy[j + 1].resize(new_s + 1);
              }
              easy[j + 1][new_s] += easy[j][s] * p_take;
              easy[j][s] *= 1.0 - p_take;
            }
          }
        }
        for (int s = 0; s < (int) dp[k].size(); s++) {
          for (int bound = dp_from[k][s]; bound < dp_from[k][s] + (int) dp[k][s].size(); bound++) {
            if (s <= bound) {
              ans += dp[k][s][bound - dp_from[k][s]];
            }
          }
        }
        for (int s = 1; s < (int) easy[k].size(); s++) {
          ans += easy[k][s];
        }
      }
      if (flag == 1) {
        reverse(rank_rounded.begin(), rank_rounded.end());
        gsea = -gsea;
      }
      memo[make_pair(k, gsea)] = make_pair(ans, total_removed);
    }
    out_f << k << "\t" << gsea << "\t" << ans << "\t" << total_removed << endl;
  }
  fclose(stdin);
  return 0;
}
