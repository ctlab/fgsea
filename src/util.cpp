#include "util.h"
#include <vector>
#include <random>

// generate k numbers from [a, b] - closed interval
// a should be non-negative, usually 0 or 1
std::vector<int> combination(const int &a, const int &b, const int &k, std::mt19937& rng) {
    // std::uniform_int_distribution<int> uni(a, b);
    uid_wrapper uni(a, b, rng);
    std::vector<int> v;
    v.reserve(k);
    std::vector<char> used(b + 1);
    for (int i = 0; i < k; i++) {
        for (int j = 0; j < 100; j++) { // average < 2
            int x = uni();
            if (!used[x]) {
                v.push_back(x);
                used[x] = true;
                break;
            }
        }
    }
    return v;
}

//*
uid_wrapper::uid_wrapper(int _from, int _to, std::mt19937& _rng) : from(_from), len(_to - _from + 1), rng(_rng) {
    unsigned maxVal = rng.max();
    completePart = maxVal - maxVal % len;
}

int uid_wrapper::operator()() {
    unsigned x;
    do {
        x = rng();
    } while (x >= completePart);

    return from + x % len;
}

/*/

uid_wrapper::uid_wrapper(int _from, int _to, std::mt19937& _rng) : rng(_rng), uid(_from, _to) {}

int uid_wrapper::operator()() {
    return uid(rng);
}
//*/