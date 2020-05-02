#include "esCalculation.h"
#include <algorithm>


double calcES(const vector<double> &ranks, const vector<int> &p, double NS) {
    // p must be sorted
    int n = (int) ranks.size();
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
        cur += q2 * ranks[pos];
        if (abs(cur) > abs(res)) {
            res = cur;
        }
        last = pos;
    }
    return res;
}

double calcES(const vector<double> &ranks, const vector<int> &p) {
    // p must be sorted
    double NS = 0.0;
    for (int pos : p) {
        NS += ranks[pos];
    }
    return calcES(ranks, p, NS);
}

double calcPositiveES(const vector<double> &ranks, const vector<int> &p, double NS) {
    // p must be sorted
    int n = (int) ranks.size();
    int k = (int) p.size();
    double res = 0.0;
    double cur = 0.0;
    double q1 = 1.0 / (n - k);
    double q2 = 1.0 / NS;
    int last = -1;
    for (int pos : p) {
        cur += q2 * ranks[pos] - q1 * (pos - last - 1);
        res = max(res, cur);
        last = pos;
    }
    return res;
}

double calcPositiveES(const vector<double> &ranks, const vector<int> &p) {
    // p must be sorted
    double NS = 0.0;
    for (int pos : p) {
        NS += ranks[pos];
    }
    return calcPositiveES(ranks, p, NS);
}


/*
bool compareStat(const vector<double> &ranks, const vector<int> &p, double NS, double bound){
    // p must be sorted
    int n = (int) ranks.size();
    int k = (int) p.size();
    double cur = 0.0;
    double q1 = 1.0 / (n - k);
    double q2 = 1.0 / NS;
    int last = -1;
    double left = 1;
    for (int pos : p) {
        cur += q2 * ranks[pos] - q1 * (pos - last - 1);
        if (cur > bound) {
            return true;
        }
        last = pos;
        left -= q2 * ranks[pos];
        if (cur + left < bound) {
            return false;
        }
    }
    return false;
}
*/

bool compareStat(const vector<double> &ranks, const vector<int> &p, double NS, double bound){
    // p must be sorted
    int n = (int) ranks.size();
    int k = (int) p.size();
    double cur = 0.0;
    double q1 = 1.0 / (n - k);
    double q2 = 1.0 / NS;
    int last = -1;
    for (int pos : p) {
        cur += q2 * ranks[pos] - q1 * (pos - last - 1);
        if (cur > bound) {
            return true;
        }
        last = pos;
    }
    return false;
}