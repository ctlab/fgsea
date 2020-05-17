#include "util.h"
#include <vector>
#include <random>

// generate k numbers from [a, b] - closed interval
// a should be non-negative, usually 0 or 1
std::vector<int> combination(const int &a, const int &b, const int &k, std::mt19937& rng) {
    std::uniform_int_distribution<int> uni(a, b);
    std::vector<int> v;
    v.reserve(k);
    std::vector<char> used(b + 1);
    for (int i = 0; i < k; i++) {
        for (int j = 0; j < 100; j++) { // average < 2
            int x = uni(rng);
            if (!used[x]) {
                v.push_back(x);
                used[x] = true;
                break;
            }
        }
    }
    return v;
}
